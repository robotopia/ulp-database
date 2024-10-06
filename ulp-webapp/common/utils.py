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
from scipy.optimize import curve_fit

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

    lc_dm_freq = lightcurve.dm_freq or np.inf # Since we stipulated that an empty field (i.e. null value) "means" infinite frequency

    total_offset = D * (lightcurve.dm - working_ephemeris.dm) / lc_dm_freq**2

    # If the above ^^^ formula is not clear, see the following "spelled out" code:
    '''
    ###################
    # Recover the time offset needed to put the lightcurve back to its "native"
    # frequency
    native_offset = D * lightcurve.dm / lc_dm_freq**2
    
    # An now the time offset to get it to infinite frequency assuming the
    # "working" DM
    inf_offset = D * working_ephemeris.dm / lc_dm_freq**2

    # We define "dm correction" as the amount you have to *add* in order to
    # get the corrected amount
    total_offset = native_offset - inf_offset
    ##################
    '''

    return total_offset


def scale_to_frequency(freq_MHz, S_freq, freq_target_MHz, alpha, q=0):
    '''
    Compare scale_flux() in common.js
    '''
    # Convert things to GHz, as that is what was assumed to be used to produce
    # alpha and q
    f = freq_MHz / 1e3
    lnf = np.log(f)
    f_target = freq_target_MHz / 1e3
    lnf_target = np.log(f_target)

    S1GHz = S_freq / (f**alpha * np.exp(q*lnf**2))
    S_target = S1GHz * f_target**alpha * np.exp(q*lnf_target**2)

    return S_target


def calc_and_create_toa(pulse_number, template, freq_target_MHz=1000):

    # Get a shorthand variable for the (w)orking (e)phemeris
    we = template.working_ephemeris

    # Get all lightcurves in this pulse number
    times, values = we.extract_lightcurves_from_pulse_number(pulse_number, freq_target_MHz)

    # Because of the awkwardness of dealing with lightcurves with generally different
    # sampling rates, we simply fit the template to the points, rather than
    # trying to use a method based on cross-correlation.
    #
    # Keep in mind that the template defines a *shape*, so any fitting function should have
    # an "amplitude" as a free parameter. This fitted amplitude can be stored with the ToA
    # as a fitted parameter

    # Convert the times to phases to prepare for fitting the template
    _, phases = we.fold(times)

    def template_func(phase, ph_offset, ampl):
        return ampl*template.values(phases - ph_offset)

    p0 = (0.0, np.max(values))
    bounds = ((-np.inf, 0.0), (np.inf, np.inf))
    popt, pcov = curve_fit(template_func, phases, values, p0=p0, bounds=bounds)

    # Unpack the fitted values
    ph_offset, ampl = popt
    ph_offset_err, ampl_err = np.sqrt(np.diag(pcov))

    # Convert the fitted phase offsets back to ToA units
    toa_mjd = we.unfold(pulse_number, ph_offset)
    toa_err_s = ph_offset_err * we.p0

    # Pack the results into a bona fide ToA
    toa = data_models.Toa(
        pulse_number=pulse_number,
        template=template,
        toa_mjd=toa_mjd,
        toa_err_s=toa_err_s,
        ampl=ampl,
        ampl_err=ampl_err,
        ampl_ref_freq=freq_target_MHz,
    )
    toa.save()

    # Return the new Toa
    return toa
