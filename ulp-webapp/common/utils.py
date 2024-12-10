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


def fit_toa(pulse, template, baseline_degree=None):

    # Get lightcurve for this pulse
    lc = pulse.lightcurve
    times = lc.bary_times(dm=0.0)
    t0 = times[0]
    values = lc.values()

    # Because of the awkwardness of dealing with lightcurves with generally different
    # sampling rates, we simply fit the template to the points, rather than
    # trying to use a method based on cross-correlation.
    #
    # Keep in mind that the template defines a *shape*, so any fitting function should have
    # an "amplitude" as a free parameter. This fitted amplitude can be stored with the ToA
    # as a fitted parameter

    # BTW, curve_fit doesn't like it when MJDs are used, since the times to be fitted are
    # very low down in precision. Better to do everything in terms of time since the start
    # and then add it all back afterwards. That's why t0 is preserved above.

    # Set up the initial values
    max_idx = np.argmax(values)
    p0 = (times[max_idx] - t0, values[max_idx])

    if baseline_degree == 0:
        def template_func(time, toa_mjd, ampl, baseline_level):
            return ampl*template.values(time - toa_mjd) + baseline_level
        p0 += (0.0,)
        bounds = ((-np.inf, 0.0, -np.inf), (np.inf, np.inf, np.inf))
    elif baseline_degree == 1:
        def template_func(time, toa_mjd, ampl, baseline_level, baseline_slope):
            return ampl*template.values(time - toa_mjd) + baseline_level + baseline_slope*(time - toa_mjd)
        p0 += (0.0, 0.0,)
        bounds = ((-np.inf, 0.0, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf))
    else:
        def template_func(time, toa_mjd, ampl):
            return ampl*template.values(time - toa_mjd)
        bounds = ((-np.inf, 0.0), (np.inf, np.inf))

    popt, pcov = curve_fit(template_func, times - t0, values, p0=p0, bounds=bounds)

    # Unpack the fitted values
    if baseline_degree == 0:
        toa_mjd, ampl, baseline_level = popt
        toa_mjd += t0
        toa_mjd_err, ampl_err, baseline_level_err = np.sqrt(np.diag(pcov))
        baseline_slope = None
    elif baseline_degree == 1:
        toa_mjd, ampl, baseline_level, baseline_slope = popt
        toa_mjd += t0
        toa_mjd_err, ampl_err, baseline_level_err, baseline_slope_err = np.sqrt(np.diag(pcov))
    else:
        toa_mjd, ampl = popt
        toa_mjd += t0
        toa_mjd_err, ampl_err = np.sqrt(np.diag(pcov))
        baseline_level = None
        baseline_slope = None

    # Convert the fitted phase offsets back to ToA units
    toa_err_s = toa_mjd_err * 86400.0

    # Pack the results into a bona fide ToA
    # Check if there is already a ToA for this pulse number
    toas = data_models.Toa.objects.filter(
        pulse=pulse,
        template=template,
    )

    if toas.exists():
        toa = toas.first()
        toa.toa_mjd = toa_mjd
        toa.err_s = toa_err_s
        toa.ampl = ampl
        toa.ampl_err = ampl_err
        toa.ampl_ref_freq = 1000.0
        toa.baseline_level = baseline_level
        toa.baseline_slope = baseline_slope
    else:
        toa = data_models.Toa(
            pulse=pulse,
            template=template,
            toa_mjd=toa_mjd,
            toa_err_s=toa_err_s,
            ampl=ampl,
            ampl_err=ampl_err,
            ampl_ref_freq=1000.0,
            baseline_level=baseline_level,
            baseline_slope=baseline_slope,
        )

    toa.save()

    # Return the new Toa
    return toa

