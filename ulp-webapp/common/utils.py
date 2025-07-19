from django.db.models import Q
from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.apps import apps
from django.contrib.auth.models import User, Group, AnonymousUser
from django.contrib.auth.decorators import login_required
import json
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import data.models as data_models
import numpy as np
from scipy.optimize import curve_fit

def permitted_to_view_filter(queryset, user):

    if user is None or user.is_anonymous:
        return queryset.none()

    return queryset.filter(
        Q(owner=user) |
        Q(published_in__isnull=False) |
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


@login_required
def update_permissions(request):

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


# Collection of fitting functions, designed to work with curve_fit, that have
# different combinations of parameters to fit for.
#

def fit_ephemeris_pepoch_p0(barycentred_and_dedispersed_toa_mjds, pepoch_mjd, p0_s):
    p0_d = p0_s/86400.0
    pulse_phases = (barycentred_and_dedispersed_toa_mjds - pepoch_mjd)/p0_d
    nearest_pulse_phases = np.round(pulse_phases)
    return nearest_pulse_phases*p0_d + pepoch_mjd
