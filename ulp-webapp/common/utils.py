from django.db.models import Q
from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.apps import apps
from django.contrib.auth.models import User, Group
import json
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import data.models as data_models
import numpy as np

def permitted_to_view_filter(queryset, user):

    return queryset.filter(
        Q(owner=user) |
        Q(can_edit_groups__user=user) |
        Q(can_edit_users=user) |
        Q(can_view_groups__user=user) |
        Q(can_view_users=user)
    ).distinct()


def permitted_to_edit_filter(queryset, user):

    return queryset.filter(
        Q(owner=user) |
        Q(can_edit_groups__user=user) |
        Q(can_edit_users=user)
    ).distinct()


def permitted_to_delete_filter(queryset, user):

    return queryset.filter(Q(owner=user)).distinct()


def update_permissions(request):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=404)

    # Turn the data into a dictionary
    try:
        data = json.loads(request.body.decode('utf-8'))
        app_name = data['app']
        model_name = data['model']
        pk = int(data['pk'])
        group_or_user = data['group_or_user']
        name = data['name']
        permission_type = data['permission_type']
        permit = bool(data['permit'])

        model = apps.get_model(app_name, model_name)
        obj = get_object_or_404(model, pk=pk)

    except:
        return HttpResponse(status=400)

    if not request.user == obj.owner:
        return HttpResponse(status=401)

    if group_or_user == 'group':
        group = get_object_or_404(Group, name=name)
        if permission_type == 'view':
            if permit:
                obj.can_view_groups.add(group)
            else:
                obj.can_view_groups.remove(group)
        elif permission_type == 'edit':
            if permit:
                obj.can_edit_groups.add(group)
            else:
                obj.can_edit_groups.remove(group)
    elif group_or_user == 'user':
        user = get_object_or_404(User, username=username)
        if permission_type == 'view':
            if permit:
                obj.can_view_users.add(user)
            else:
                obj.can_view_users.remove(user)
        elif permission_type == 'edit':
            if permit:
                obj.can_edit_users.add(user)
            else:
                obj.can_edit_users.remove(user)

    return HttpResponse(status=200)


def barycentre(ulp, times, location):
    '''
    ulp is a Ulp object
    times is an array-like object of MJDs (can be Time array)
    '''

    ra, dec = [
        data_models.EphemerisMeasurement.objects.filter(
            ephemeris_parameter__tempo_name=param,
            measurement__ulp=ulp
            ).first().measurement.quantity for param in ['RAJ', 'DECJ']
    ]

    direction = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    location_times = Time(times, format='mjd', scale='utc', location=location) # Tag times with location
    bc_correction = location_times.light_travel_time(direction, ephemeris='jpl') # Calculate correction
    bc_times = (location_times.tdb + bc_correction).value # Apply correction, and convert back to MJD values
    return bc_times


# Compare fold() in common.js
def fold(mjds, working_ephemeris):
    pepoch = working_ephemeris.pepoch
    period = working_ephemeris.p0

    pulses_phases = (mjds - pepoch) / (period/86400.0)
    phases, pulses = np.modf(pulses_phases + 0.5)
    phases -= 0.5

    return pulses, phases


def unfold(pulses, phases, working_ephemeris):
    pepoch = working_ephemeris.pepoch
    period = working_ephemeris.p0
    mjds = pepoch + period*(pulses + phases)
    return mjds


def pulses_from_pulse_number(pulse_number, working_ephemeris):
    mjd_start, mjd_end = unfold(pulse_number, np.array([-0.5, 0.5]), working_ephemeris)
    return models.Pulse.objects.filter(mjd_start__gte=mjd_start, mjd_start__lte=mjd_end)


def dm_correction(lightcurve, working_ephemeris):
    '''
    This functions solves the problem of when the lightcurve was dedispersed using
    some particular DM and some arbitrary reference frequency, whereas the working
    ephemeris is assuming a different DM. We can't change the channel dispersion
    delays, but we can find the absolute time offset to approximate how (and when)
    the lightcurve appears if it had been dedispersed to infinite frequency.

    Returns time offset in units of days, in order that it can be easily to applied
    to ToAs, which are in MJD.
    '''

    D = 4.148808e3/86400 # Dispersion constant in the appropriate units

    total_offset = D * (lightcurve.dm - working_ephemeris.dm) / lightcurve.dm_freq**2

    # If the above ^^^ formula is not clear, see the following "spelled out" code:
    '''
    ###################
    # Recover the time offset needed to put the lightcurve back to its "native"
    # frequency
    native_offset = D * lightcurve.dm / lightcurve.dm_freq**2
    
    # An now the time offset to get it to infinite frequency assuming the
    # "working" DM
    inf_offset = D * working_ephemeris.dm / lightcurve.dm_freq**2

    # We define "dm correction" as the amount you have to *add* in order to
    # get the corrected amount
    total_offset = native_offset - inf_offset
    ##################
    '''

    return total_offset


def calc_and_create_toa(pulse_number, template, working_ephemeris):

    # Get the pulses we'll be working with
    pulses = pulses_from_pulse_number(pulse_number, working_ephemeris)

    # Extract the lightcurves and shift them to infinite frequency
    times = np.concatenate([p.lightcurve.bary_times() + dm_correction(lightcurve, working_ephemeris) for p in pulses])
    values = np.concatenate([p.lightcurve.values() for p in values])

    # If there are more than one pulse, then there is more work to be done in getting
    # a straightforward array which we can use for cross correlating with the template.
    # Current plan is:
    # 1. Order the samples
    # 2. Interpolate and resample (scipy) to some fixed sampling time
    # For now, however, this will remain a job for tomorrow me.
    # TODO
    if len(pulses) > 1:
        raise ValueError("Currently do not support creating ToAs when multiple 'pulses' are present in the same pulse")

    # Convert the times into phases
    _, phases = fold(times, working_ephemeris)

    # Create a period-sized template at the same sampling rate as the data
    # (TODO: this needs to be more intelligent when multiple pulses are included)
    dph = pulses[0].lightcurve.dt/working_ephemeris.p0 # (dph = "Delta phase")

    # Construct a template that ....
    # TODO Finish me!
