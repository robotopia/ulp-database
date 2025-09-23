from django.shortcuts import render, get_object_or_404, redirect
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.db.models import Q, Value, BooleanField
from django.urls import reverse
from django.core.exceptions import ValidationError
from . import models
from . import serializers
from common.utils import *
from django.utils.timezone import now
from urllib.parse import urlencode

from rest_framework.authentication import TokenAuthentication
from rest_framework.permissions import IsAuthenticated
from rest_framework.decorators import api_view, authentication_classes, permission_classes
import rest_framework.exceptions as rest_exceptions

from published import models as published_models
from published.views import get_accessible_measurements

import numpy as np
from scipy.optimize import curve_fit
import json
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz, get_sun
from astropy.constants import c
from astropy.time import Time
from astropy.table import QTable
from decimal import Decimal
import io
import re

from spiceypy.spiceypy import spkezr, furnsh, j2000, spd, unload
import pandas as pd
from plotly.offline import plot
import plotly.express as px
#import plotly.graph_objects as go
#from plotly.subplots import make_subplots

# Curate list of supported telescope sites
EarthLocation.__hash__ = lambda loc: hash((loc.x, loc.y, loc.y))
sites_reversed = {EarthLocation.of_site(site_name): site_name.lower() for site_name in reversed(EarthLocation.get_site_names())}
sites = dict(zip(sites_reversed.values(), sites_reversed.keys()))

# Add extra sites not in AstroPy's original list.
# TODO -------------------  get better coordinates (these were just from map, e.g. Google)
sites['pushchino observatory'] = EarthLocation.from_geodetic(37.619605*u.deg, 54.825089*u.deg)
sites['atca'] = EarthLocation.from_geodetic(149.5623193*u.deg, -30.3132765*u.deg)
# END TODO ----------------

site_names = sorted(sites.keys(), key=lambda x: x.lower())


# Defined at the database level
toa_freq_units = "MHz"
toa_mjd_err_units = "d"
observation_freq_units = "MHz"
observation_duration_units = "s"
observation_start_mjd_format = "mjd"

def calc_dmdelay(dm, flo, fhi):
    return 4.148808e3 * u.s * (dm.to('pc/cm3').value) * (1/flo.to('MHz').value**2 - 1/fhi.to('MHz').value**2)

def ephemeris_to_skycoord(ephemeris):
    '''
    ephemeris_measurements_qs is a QuerySet. If it contains both RA and DEC,
    a SkyCoord object will be constructed therefrom.
    '''
    try:
        coord = SkyCoord(ra=ephemeris['ra'], dec=ephemeris['dec'], unit=(u.deg, u.deg), frame='icrs')
    except:
        coord = None

    return coord


def calc_pulse_phase(time, ephemeris):
    '''
    time is an astropy Time object
    ephemeris['pepoch'] should be in days
    ephemeris['p0'] should be in seconds
    ephemeris['p1'] should be dimensionless
    '''
    try:
        pepoch = Time(ephemeris['pepoch'], format='mjd')
        p0 = ephemeris['p0']*u.s
        p1 = u.Quantity(ephemeris['p1'])
        t = time - pepoch
        #return (t/p0).decompose()  # This is only good when p1 = 0
        return (t/p0 - 0.5*p1*(t/p0)**2).decompose()
    except:
        return 0.0


def calc_mjd(pulse_phase, ephemeris):
    '''
    This is the inverse of calc_pulse_phase()
    '''
    pepoch = Time(ephemeris['pepoch'], format='mjd')
    p0 = ephemeris['p0']*u.s
    p1 = u.Quantity(ephemeris['p1'])

    #return pepoch + pulse_phase*p0  # This is only good when p1 = 0
    # The "other" quadradtic formula, valid for small or zero quadratic term (i.e. p1=0)
    return 2*pulse_phase*p0 / (1 + np.sqrt(1 - 2*p1*pulse_phase)) + pepoch


def generate_toas(time_start, time_end, ephemeris):

    pulse_phase_start = calc_pulse_phase(time_start, ephemeris)
    pulse_phase_end = calc_pulse_phase(time_end, ephemeris)

    pulse_phases = np.arange(np.ceil(pulse_phase_start), pulse_phase_end)
    mjds = calc_mjd(pulse_phases, ephemeris)

    return mjds


def toa_data(user, we):

    # Get the ToAs that this user is allowed to view
    toas = permitted_to_view_filter(models.TimeOfArrival.objects.filter(ulp=we.ulp, raw_mjd__isnull=False, freq__isnull=False), user)

    # Barycentre
    mjds = Time([Time(float(toa.raw_mjd), scale='utc', format='mjd', location=sites[toa.telescope_name.lower()]) for toa in toas])
    bc_corrections = mjds.light_travel_time(we.coord, ephemeris='jpl').to('day').value

    toas_data = [
        {
            'mjd': float(toas[i].raw_mjd),
            'mjd_err': float(toas[i].mjd_err),
            'freq_MHz': float(toas[i].freq),
            'bc_correction': bc_corrections[i],
            'detail_link': reverse('toa_detail_view', args=[toas[i].pk]),
        } for i in range(len(toas))
    ]

    return toas_data


@login_required
def toa_json(request, we_pk):

    # Retrieve the selected ULP and working ephemeris
    we = get_object_or_404(models.WorkingEphemeris, pk=we_pk)

    return JsonResponse(toa_data(request.user, we), safe=False)


@login_required
def obs_data(request, we_pk):

    # Retrieve the selected ULP and working ephemeris
    we = get_object_or_404(models.WorkingEphemeris, pk=we_pk)

    # Get the ToAs that this user is allowed to view
    obss = permitted_to_view_filter(models.Observation.objects.filter(ulps=we.ulp, freq__isnull=False), request.user)

    if obss.exists():

        # Barycentre
        mjds = Time([Time(float(obs.start_mjd), scale='utc', format='mjd') + obs.duration*u.s/2 for obs in obss])
        bc_corrections = [mjds[i].light_travel_time(we.coord, ephemeris='jpl', location=sites[obss[i].telescope_name.lower()]).to('day').value for i in range(len(mjds))]

        obss_json = [
            {
                'obsid': obss[i].obsid,
                'mjd': float(mjds[i].mjd),
                'mjd_err': obss[i].duration/86400/2,
                'freq_MHz': float(obss[i].freq),
                'bc_correction': bc_corrections[i],
                'detail_link': obss[i].obsid,  # TODO: Change me!
            } for i in range(len(obss))
        ]

        return JsonResponse(obss_json, safe=False)

    # If reach here, no observations for this object
    return JsonResponse([], safe=False)


@login_required
def timing_choose_ulp_view(request):

    # Get ULPs to which they have access
    ulps = published_models.Ulp.objects.filter(
        Q(data_access_groups__user=request.user) |
        Q(whitelist_users=request.user)
    ).distinct()

    context = {
        'ulps': ulps,
        'error': request.session.pop('error', None),
    }

    return render(request, 'data/timing_choose_ulp.html', context)


def get_toa_predictions(start, end, freq, pepoch, p0, p1, dm, telescope, coord, output_toa_format='mjd', min_el=0, max_sun_el=90, p_aw=None, t0_aw=None, duration_aw=None):

    if start >= end:
        return []

    location = sites[telescope.lower()]

    # First, assume the given mjd_start and mjd_end are in fact topocentric dispersed MJDs,
    # so to get the right range, convert them to dedispersed, barycentric
    time_range = Time([start, end])
    dmdelay = calc_dmdelay(dm, freq, np.inf*u.MHz)
    time_range -= dmdelay
    time_range += time_range.light_travel_time(coord, ephemeris='jpl', location=sites[telescope.lower()])

    ephemeris = {
        'ra': coord.ra.deg,
        'dec': coord.dec.deg,
        'pepoch': pepoch,
        'p0': p0,
        'p1': p1,
        'dm': dm.to('pc cm-3').value,
    }

    predicted_barycentric_toas = generate_toas(time_range[0], time_range[1], ephemeris)
    predicted_topocentric_toas = predicted_barycentric_toas - predicted_barycentric_toas.light_travel_time(coord, ephemeris='jpl', location=sites[telescope.lower()])
    predicted_dispersed_toas = predicted_topocentric_toas + dmdelay

    altaz = coord.transform_to(AltAz(obstime=predicted_dispersed_toas, location=sites[telescope.lower()]))
    sun_altaz = [get_sun(predicted_dispersed_toas[i]).transform_to(AltAz(obstime=predicted_dispersed_toas[i], location=sites[telescope.lower()])) for i in range(len(predicted_barycentric_toas))]

    # If all necessary activity window values are provided, only choose toas
    # within the activity window
    if p_aw and t0_aw and duration_aw:
        # Apply units
        p_aw *= u.h
        t0_aw = Time(t0_aw, scale='utc', format='mjd')
        duration_aw *= u.h

        aw_number, aw_phase = np.divmod(((predicted_barycentric_toas - t0_aw)/p_aw).decompose() + 0.5, 1)
        aw_phase -= 0.5
        aw_phase_max = (0.5*duration_aw/p_aw).decompose()
        aw_phase_min = -aw_phase_max
    else:
        aw_phase = np.zeros(predicted_barycentric_toas.shape)
        aw_phase_min, aw_phase_max = -0.5, 0.5

    # Set to the requested format
    predicted_barycentric_toas.format = output_toa_format
    predicted_topocentric_toas.format = output_toa_format
    predicted_dispersed_toas.format = output_toa_format

    predicted_toas = [
        {
            'bary': predicted_barycentric_toas[i],
            'topo': predicted_topocentric_toas[i],
            'topo_disp': predicted_dispersed_toas[i],
            'elevation': int(np.round(altaz[i].alt.value)),
            'sun_elevation': int(np.round(sun_altaz[i].alt.value)),
            'aw_phase': aw_phase[i],
        } for i in range(len(predicted_barycentric_toas)) if (
            altaz[i].alt.value >= min_el and
            sun_altaz[i].alt.value < max_sun_el and
            aw_phase[i] >= aw_phase_min and
            aw_phase[i] <= aw_phase_max)
    ]

    return predicted_toas


@api_view(['GET'])
@authentication_classes([TokenAuthentication])
@permission_classes([IsAuthenticated])
def get_toa_predictions_json(request):

    ulp_name = request.GET.get('lpt')
    ulp = published_models.Ulp.objects.filter(name=ulp_name).first()
    if ulp is None:
        return JsonResponse({'message': f"Couldn't identify an LPT with the name {ulp_name}"}, status=400)

    default_date_format = 'mjd'

    input_date_format = request.GET.get('input_date_format', default_date_format)
    if input_date_format not in Time.FORMATS.keys():
        return JsonResponse({'message': f"{input_date_format} is not a valid AstroPy date format. Please select one of the following (default='{default_date_format}'):\n{', '.join(list(Time.FORMATS.keys()))}"}, status=400)

    output_toa_format = request.GET.get('output_toa_format', default_date_format)

    start = Time(request.GET.get("start"), scale='utc', format=input_date_format)
    end = Time(request.GET.get("end"), scale='utc', format=input_date_format)

    try:
        freq = float(request.GET.get("freq_MHz")) * u.MHz
    except:
        return JsonResponse({'message': f"Couldn't parse frequency '{request.GET.get('freq_MHz')}' as a valid float"}, status=400)

    telescope = request.GET.get("telescope")
    if telescope is None:
        return JsonResponse({'message': f"Must provide a telescope"}, status=400)
    if telescope not in site_names:
        return JsonResponse({'message': f"{telescope} not in AstroPy's list of 'site names'"}, status=400)

    try:
        min_el = float(request.GET.get("min_el", 0))
    except:
        return JsonResponse({'message': f"Couldn't parse min_el = '{request.GET.get('min_el')}' as a valid float"}, status=400)

    try:
        max_sun_el = float(request.GET.get("max_sun_el", 90))
    except:
        return JsonResponse({'message': f"Couldn't parse max_sun_el = '{request.GET.get('max_sun_el')}' as a valid float"}, status=400)

    # Get the working ephemerides that are available to this user
    working_ephemeris = models.WorkingEphemeris.objects.filter(ulp=ulp, owner=request.user).first()

    if not working_ephemeris:
        return JsonResponse({'message': f"No ephemerides available for {ulp.name} owned by {request.user}. Please go onto the Timing page on the website to create one."}, status=400)

    coord = working_ephemeris.coord
    pepoch = working_ephemeris.pepoch
    p0 = working_ephemeris.p0
    p1 = working_ephemeris.p1 or 0.0
    dm = working_ephemeris.dm * u.pc / u.cm**3

    predictions = get_toa_predictions(
        start,
        end,
        freq,
        pepoch,
        p0,
        p1,
        dm,
        telescope,
        coord,
        output_toa_format=output_toa_format,
        min_el=min_el,
        max_sun_el=max_sun_el,
    )
    
    for prediction in predictions:
        prediction['bary'] = prediction['bary'].to_value(output_toa_format)
        prediction['topo'] = prediction['topo'].to_value(output_toa_format)
        prediction['topo_disp'] = prediction['topo_disp'].to_value(output_toa_format)

    return JsonResponse(predictions, safe=False, status=200)


@login_required
def timing_residual_view(request, pk):

    # Retrieve the selected ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    context = {'ulp': ulp}

    # Make sure the user has the permissions to view this ULP

    # Second, they have to belong to a group that has been granted access to
    # this ULP's data
    if not ulp.data_access_groups.filter(user=request.user).exists() and not request.user in ulp.whitelist_users.all():
        return redirect('timing_choose_ulp')

    # Otherwise, grant them access, and get the TOAs!
    toas = models.TimeOfArrival.objects.filter(ulp=ulp)

    # Now limit only to those that this user has view access to
    toas = permitted_to_view_filter(toas, request.user)

    if not toas.exists():
        request.session['error'] = f"No ToAs are available for {ulp}, or you do not have permission to see them."
        return redirect('timing_choose_ulp')

    # Get the working ephemerides that are available to this user
    working_ephemerides = permitted_to_view_filter(models.WorkingEphemeris.objects.filter(ulp=ulp), request.user).order_by('owner')

    # If there aren't any, make a default one for the user!
    if not working_ephemerides.filter(owner=request.user).exists():
        we = models.WorkingEphemeris.objects.filter(owner__isnull=True, ulp=ulp).first()
        if we is None:
            # No default, just create an empty one
            we = models.WorkingEphemeris(owner=request.user, ulp=ulp)
        else:
            we.pk = None # When save(), will now create a new object
            we.owner = request.user
        we.save()
        working_ephemerides = permitted_to_view_filter(models.WorkingEphemeris.objects.filter(ulp=ulp), request.user)

    context['working_ephemerides'] = working_ephemerides

    # Choose a working ephemeris to be 'selected'
    try:
        selected_working_ephemeris = working_ephemerides.get(pk=request.POST.get('working_ephemeris_pk'))
    except:
        selected_working_ephemeris = working_ephemerides.get(owner=request.user)
    if selected_working_ephemeris.covariance is None:
        selected_working_ephemeris.covariance = models.WorkingEphemerisCovariance()
        selected_working_ephemeris.covariance.save()
        selected_working_ephemeris.save()
    context['selected_working_ephemeris'] = selected_working_ephemeris

    # Pool together lists of published values to offer as options to the user
    context['published_ephemeris_values'] = {
        'ra': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Right ascension", article__isnull=False),
        'dec': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Declination", article__isnull=False),
        'pepoch': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="PEPOCH", article__isnull=False),
        'p0': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Period", article__isnull=False),
        'p1': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Period derivative", article__isnull=False),
        'dm': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Dispersion measure", article__isnull=False),
        'p_aw': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Activity window: period", article__isnull=False),
        't0_aw': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Activity window: reference epoch", article__isnull=False),
        'duration_aw': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Activity window: duration", article__isnull=False),
    }

    # Construct a dictionary out of the ephemeris
    ephemeris = {
        'ra': selected_working_ephemeris.ra,
        'dec': selected_working_ephemeris.dec,
        'pepoch': selected_working_ephemeris.pepoch,
        'p0': selected_working_ephemeris.p0,
        'p1': selected_working_ephemeris.p1 or 0.0,
        'dm': selected_working_ephemeris.dm,
        'p_aw': selected_working_ephemeris.p_aw or 0.0,
        't0_aw': selected_working_ephemeris.t0_aw or 0.0,
        'duration_aw': selected_working_ephemeris.duration_aw or 0.0,
    }
    coord = selected_working_ephemeris.coord

    # Do something similar for toas that are not yet dedispersed, but for which a frequency is given
    if 'dm' in ephemeris.keys() and ephemeris['dm'] is not None:
        dm = ephemeris['dm'] * u.pc / u.cm**3
        for toa in toas:
            if toa.freq is not None: # and toa.dedispersed == False:
                delay = calc_dmdelay(dm, toa.freq*u.MHz, np.inf*u.MHz)
                toa.mjd -= Decimal(delay.to('d').value)
                #toa.dedispersed = True
                #toa.save()

    output_toa_format = 'mjd'
    min_el = 0.0
    max_sun_el = 90.0
    mjd_start_format = 'iso'
    mjd_end_format = 'iso'

    if request.method == "POST":

        # Get form values
        pepoch = float(request.POST.get('pepoch', '0.0'))
        p0 = float(request.POST.get('p0', '0.0'))
        p1 = float(request.POST.get('p1', '0.0'))
        dm = float(request.POST.get('dm', '0.0'))
        ra = request.POST.get('ra', '00:00:00')
        dec = request.POST.get('dec', '00:00:00')
        p_aw = float(request.POST.get('p_aw', '0.0'))
        t0_aw = float(request.POST.get('t0_aw', '0.0'))
        duration_aw = float(request.POST.get('duration_aw', '0.0'))

        mjd_start_format = request.POST.get('mjd-start-format')
        mjd_end_format = request.POST.get('mjd-end-format')
        mjd_start = Time(request.POST.get('mjd-start'), format=mjd_start_format)
        mjd_end = Time(request.POST.get('mjd-end'), format=mjd_end_format)
        mjd_range = Time([mjd_start, mjd_end])

        context['mjd_start'] = mjd_start.to_value(mjd_start_format)
        context['mjd_end'] = mjd_end.to_value(mjd_start_format)

        telescope = request.POST.get('telescope')
        try:
            min_el = float(request.POST.get('minimum-elevation'))
        except:
            min_el = 0.0
        try:
            max_sun_el = float(request.POST.get('maximum-sun-elevation'))
        except:
            max_sun_el = 90.0

        context['selected_telescope'] = telescope

        mjd_dispersion_frequency = float(request.POST.get('mjd-dispersion-frequency'))
        output_toa_format = request.POST.get('output-toa-format')

        # Populate the ephemeris from the form values
        ephemeris['pepoch'] = pepoch
        ephemeris['p0'] = p0
        ephemeris['p1'] = p1
        ephemeris['dm'] = dm
        ephemeris['ra'] = ra
        ephemeris['dec'] = dec
        ephemeris['p_aw'] = p_aw
        ephemeris['t0_aw'] = t0_aw
        ephemeris['duration_aw'] = duration_aw

        # If they've also provided other form values, make a table of predicted values
        if mjd_start is not None and mjd_end is not None and mjd_dispersion_frequency is not None and pepoch is not None and p0 is not None and telescope is not None:
            context['predicted_toas'] = get_toa_predictions(
                mjd_range[0],                        # Start of time range
                mjd_range[-1],                       # End of time range
                mjd_dispersion_frequency*u.MHz,      # Observing frequency
                pepoch,                              # PEPOCH (MJD)
                p0,                                  # Period
                p1,                                  # Period derivative
                ephemeris['dm']*u.pc/u.cm**3,        # DM
                telescope,                           # Telescope string
                coord,                               # Target source coordinates
                output_toa_format=output_toa_format, # Desired output ToA format
                min_el=min_el,                       # Minimum source elevation
                max_sun_el=max_sun_el,               # Maximum Sun elevation
                p_aw=p_aw,                           # Activity window period (hours)
                t0_aw=t0_aw,                         # Activity window reference epoch (MJD)
                duration_aw=duration_aw,             # Activity window duration (hours)
            )

        context['mjd_dispersion_frequency'] = mjd_dispersion_frequency

    context['ephemeris'] = ephemeris
    context['time_formats'] = list(Time.FORMATS.keys())
    context['time_formats'].sort()

    # Get available published periods
    periods = published_models.Measurement.objects.filter(
        parameter__name="Period", # Hard code this specific parameter name
        ulp=ulp,
    )
    periods = periods.filter(
        Q(article__isnull=False) |  # It's published, and therefore automatically accessible by everyone
        Q(owner=request.user) |  # The owner can always see their own measurements
        Q(ulp__whitelist_users=request.user) |  # Allow whitelisted users
        Q(access=published_models.Measurement.ACCESS_PUBLIC) |  # Include measurements explicitly marked as public
        (Q(access=published_models.Measurement.ACCESS_GROUP) &  # But if it's marked as group-accessible...
         Q(access_groups__in=request.user.groups.all()))  # ...then the user must be in of the allowed groups.
    )
    context['periods'] = periods

    # Calculate some sensible initial plot dimensions
    mjds = Time([float(toa.mjd) for toa in toas], format='mjd')
    xdata_min = np.min(mjds).mjd
    xdata_max = np.max(mjds).mjd
    xdata_range = xdata_max - xdata_min

    plot_specs = {
        'xmin': xdata_min - 0.05*xdata_range,  # With an extra margin buffer
        'xmax': xdata_max + 0.05*xdata_range,
        'xrange': xdata_range,
    }
    context['plot_specs'] = plot_specs

    # Generate list of output TOA formats:
    context['selected_output_toa_format'] = output_toa_format
    context['telescopes'] = site_names
    context['min_el'] = min_el
    context['max_sun_el'] = max_sun_el
    context['mjd_start_format'] = mjd_start_format
    context['mjd_end_format'] = mjd_start_format

    context['covariance'] = serializers.WorkingEphemerisCovarianceSerializer().serialize(
        models.WorkingEphemerisCovariance.objects.filter(pk=selected_working_ephemeris.covariance.pk)
    )

    return render(request, 'data/timing_residuals.html', context)


@login_required
def toa_detail_view(request, pk):

    # Retrieve the selected ULP
    toa = get_object_or_404(models.TimeOfArrival, pk=pk)

    context = {
        'toa': toa,
        'freq_units': toa_freq_units,
    }

    return render(request, 'data/toa_detail.html', context)

@login_required
def toas_view(request, pk):

    # Retrieve the selected ToAs from the specified ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)
    toas = ulp.times_of_arrival.all()

    # Now limit only to those that this user has view access to
    toas = permitted_to_view_filter(toas, request.user)

    # Make sure the requested display units are dimensionally correct; if not, use default
    freq_units = "MHz"
    mjd_err_units = "us"

    # Annotate ToAs according to whether the user can edit it or not
    # See, e.g., https://stackoverflow.com/questions/41354910/how-to-annotate-the-result-of-a-model-method-to-a-django-queryset
    for toa in toas:
        toa.editable = toa.can_edit(request.user)
        # Also switch units to whatever is requested
        toa.freq = (toa.freq * u.Unit(toa_freq_units)).to(freq_units).value if toa.freq else None
        toa.bw = (toa.bw * u.Unit(toa_freq_units)).to(freq_units).value if toa.bw else None
        toa.mjd_err = (toa.mjd_err * u.Unit(toa_mjd_err_units)).to(mjd_err_units).value

    context = {
        "ulp": ulp,
        "toas": toas,
        "telescopes": site_names,
        "freq_units": freq_units,
        "mjd_err_units": mjd_err_units,
    }

    return render(request, 'data/toas.html', context)


@login_required
def observations_view(request, pk):

    # Retrieve the selected Observations from the specified ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)
    observations = ulp.observations.all()

    # Now limit only to those that this user has view access to
    observations = permitted_to_view_filter(observations, request.user)

    # Apply filters
    try:
        from_date = Time(request.GET.get('from_date'), scale='utc', format=request.GET.get('time_format'))
    except:
        from_date = None
    if from_date:
        observations = observations.filter(start_mjd__gte=from_date.mjd)

    try:
        to_date = Time(request.GET.get('to_date'), scale='utc', format=request.GET.get('time_format'))
    except:
        to_date = None
    if to_date:
        observations = observations.filter(start_mjd__lte=to_date.mjd)

    # Paging!
    try:
        page = int(request.GET.get("page"))
    except:
        page = 1

    try:
        page_size = int(request.GET.get("page_size"))
    except:
        page_size = 50

    start_idx = (page - 1)*page_size
    end_idx = start_idx + page_size
    npages = len(observations) // page_size + 1

    if start_idx >= len(observations) or start_idx < 0:
        observations = models.Observation.objects.none()
    elif end_idx > len(observations):
        observations = observations[start_idx:]
    else:
        observations = observations[start_idx:end_idx]

    # Make sure the requested display units are dimensionally correct; if not, use default
    freq_units = "MHz"
    bw_units = "MHz"
    duration_units = "s"
    start_mjd_format = "mjd"

    # Annotate observations according to whether the user can edit it or not
    # See, e.g., https://stackoverflow.com/questions/41354910/how-to-annotate-the-result-of-a-model-method-to-a-django-queryset
    for observation in observations:
        observation.editable = observation.can_edit(request.user)
        # Also switch units to whatever is requested
        observation.freq = (observation.freq * u.Unit(observation_freq_units)).to(freq_units).value if observation.freq else None
        observation.bw = (observation.bw * u.Unit(observation_freq_units)).to(bw_units).value if observation.bw else None
        observation.duration = (observation.duration * u.Unit(observation_duration_units)).to(duration_units).value if observation.duration else None

    context = {
        "ulp": ulp,
        "observations": observations,
        "telescopes": site_names,
        "column_formats": {
            "duration_units": duration_units,
            "freq_units": freq_units,
            "bw_units": bw_units,
            "start_mjd_format": start_mjd_format,
        },
        "pages": {
            "first": 1,
            "prev": page-1 if page > 1 else 1,
            "this": page,
            "next": page+1 if page < npages else npages,
            "last": npages,
            "size": page_size,
        },
        "filters": {
            'time_formats': sorted(list(Time.FORMATS.keys())),
            'selected_time_format': request.GET.get('time_format', 'mjd'),
            'from_date': request.GET.get('from_date'),
            'to_date': request.GET.get('to_date'),
        }
    }

    return render(request, 'data/observations.html', context)


def add_or_update_pulse(request, pk):

    lightcurve = get_object_or_404(models.Lightcurve, pk=pk)

    if request.method == "POST":

        # I should be validating data here...

        # Retrieve the peak value in the given range
        def get_peak(lc, mjd_start, mjd_end):
            sample_number_start = (mjd_start - lc.t0)/(lc.dt/86400)
            sample_number_end = (mjd_end - lc.t0)/(lc.dt/86400)
            points = lc.points.filter(
                sample_number__gte=sample_number_start,
                sample_number__lte=sample_number_end,
            )
            peak_value_idx = int(np.argmax([p.value for p in points]))
            peak_value_Jy = points[peak_value_idx].value
            peak_value_sample_number = points[peak_value_idx].sample_number
            peak_value_MJD = peak_value_sample_number*(lc.dt/86400) + lc.t0

            return peak_value_MJD, peak_value_Jy

        if 'pulse_id' in request.POST.keys():
            # We're updating an existing pulse
            pulse = get_object_or_404(models.Pulse, pk=int(request.POST['pulse_id']))

            peak_value_MJD, peak_value_Jy = get_peak(
                lightcurve,
                float(request.POST['mjd_start']),
                float(request.POST['mjd_end']),
            )

            if request.POST['action'] == "Save":
                pulse.mjd_start = request.POST['mjd_start']
                pulse.mjd_end = request.POST['mjd_end']
                pulse.peak_value_Jy = peak_value_Jy
                pulse.peak_value_MJD = peak_value_MJD
                pulse.tags = request.POST['tags']
                pulse.save()
            elif request.POST['action'] == "Delete":
                print(f"about to delete {pulse}...")
                pulse.delete()

        else:
            # We're creating a new pulse
            peak_value_MJD, peak_value_Jy = get_peak(
                lightcurve,
                float(request.POST['mjd_start']),
                float(request.POST['mjd_end']),
            )

            pulse = models.Pulse(
                lightcurve=lightcurve,
                mjd_start=request.POST['mjd_start'],
                mjd_end=request.POST['mjd_end'],
                peak_value_Jy=peak_value_Jy,
                peak_value_MJD=peak_value_MJD,
                tags=request.POST['tags'],
            )

            pulse.save()

    return redirect('lightcurve_view', pk=pk)


def convert_units(request):

    # Turn the data into a dictionary
    data = json.loads(request.body.decode('utf-8'))

    # Get the values to be converted, and the units conversions to be effected
    try:
        values = data['values']
        from_unit = data['from_unit']
        to_unit = data['to_unit']

        # Do the conversion
        # astropy.units and/or numpy can't handle None's being in the array, so convert to NaNs first
        values = [value if value != None else np.nan for value in values]
        quantities = (values * u.Unit(from_unit)).to(to_unit)
        new_values = [value if np.isfinite(value) else None for value in quantities.value]

    except Exception as err:
        return HttpResponse(str(err), status=400)

    return JsonResponse(new_values, safe=False, status=200)


def convert_time_format(request):

    # Turn the data into a dictionary
    data = json.loads(request.body.decode('utf-8'))

    # Get the values to be converted, and the units conversions to be effected
    try:
        values = data['values']
        from_format = data['from_format']
        to_format = data['to_format']

        # Do the conversion
        # astropy.Time can't handle None's being in the array, so have to do this one value at a time
        new_values = [getattr(Time(value, scale='utc', format=from_format), to_format) if value != None else None for value in values]

    except Exception as err:
        return HttpResponse(str(err), status=400)

    return JsonResponse(new_values, safe=False, status=200)


@login_required
def act_on_toas(request, pk):

    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Parse out the ToA primary keys from the selected checkbox names
    toa_pks = [int(key[3:]) for key in request.POST.keys() if key.startswith('cb_')]

    if request.POST.get('action_on_selected') == 'delete':

        # Get only the ToAs the user has permission to dete
        toas = permitted_to_delete_filter(models.TimeOfArrival.objects.filter(pk__in=toa_pks), request.user)

        toas.delete()

    return redirect('toas_view', pk=pk)


@login_required
def act_on_observations(request, pk):

    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Parse out the ToA primary keys from the selected checkbox names
    observation_pks = [int(key[3:]) for key in request.POST.keys() if key.startswith('cb_')]

    if request.POST.get('action_on_selected') == 'delete':

        # Get only the observations the user has permission to dete
        observations = permitted_to_delete_filter(models.Observation.objects.filter(pk__in=observation_pks), request.user)

        observations.delete()

    return redirect('observations_view', pk=pk)


@login_required
def add_toa(request, pk):

    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    try:
        freq = float(request.POST.get('freq'))
        raw_mjd = float(request.POST.get('raw_mjd'))
        mjd_err = float(request.POST.get('mjd_err'))
        telescope_name = request.POST.get('telescope_name')
        freq_units = request.POST.get('freq_units')
        mjd_err_units = request.POST.get('mjd_err_units')

        freq *= u.Unit(freq_units)
        mjd_err *= u.Unit(mjd_err_units)
    except Exception as err:
        return HttpResponse(str(err), status=400)

    toa = models.TimeOfArrival(
        owner=request.user,
        ulp=ulp,
        freq=freq.to(toa_freq_units).value,
        raw_mjd=raw_mjd,
        mjd_err=mjd_err.to(toa_mjd_err_units).value,
        telescope_name=telescope_name,
        mjd=raw_mjd,
        barycentred=False,
        dedispersed=False,
    )

    toa.save()

    return redirect('toas_view', pk=pk)


@login_required
def add_observation(request, pk):

    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    try:
        telescope_name = request.POST.get('telescope_name')
        freq = float(request.POST.get('freq'))
        bw = float(request.POST.get('bw'))
        start_mjd = request.POST.get('start_mjd')
        duration = float(request.POST.get('duration'))
        freq_units = request.POST.get('freq_units')
        bw_units = request.POST.get('bw_units')
        duration_units = request.POST.get('duration_units')
        start_mjd_format = request.POST.get('start_mjd_format')

        freq *= u.Unit(freq_units)
        bw *= u.Unit(bw_units)
        duration *= u.Unit(duration_units)
        start_mjd = Time(start_mjd, scale='utc', format=start_mjd_format).mjd

    except Exception as err:
        return HttpResponse(str(err), status=400)

    observation = models.Observation(
        owner=request.user,
        telescope_name=telescope_name,
        freq=freq.to(observation_freq_units).value,
        bw=bw.to(observation_freq_units).value,
        start_mjd=start_mjd,
        duration=duration.to(observation_duration_units).value,
    )

    observation.save()
    observation.ulps.add(ulp)

    return redirect('observations_view', pk=pk)


@login_required
def update_toa(request):

    # Turn the data into a dictionary
    data = json.loads(request.body.decode('utf-8'))

    # Get the relevant TOA object
    toa = get_object_or_404(models.TimeOfArrival, pk=data['pk'])

    # Second, they have to belong to a group that has been granted edit privileges
    if not toa.can_edit(request.user):
        return HttpResponse(status=403)

    # Set the field value
    value = data['value']
    unit = u.Unit(data['unit'])
    if data['field'] in ["mjd_err"]:
        value = (float(value)*unit).to(toa_mjd_err_units).value
    if data['field'] in ["freq", "bw"]:
        value = (float(value)*unit).to(toa_freq_units).value
    setattr(toa, data['field'], value)

    # Save the result to the database
    try:
        toa.save()
    except ValidationError as err:
        return HttpResponse(str(err), status=400)

    # Return "all is well"
    return HttpResponse(status=200)


@login_required
def update_observation(request):

    # Turn the data into a dictionary
    data = json.loads(request.body.decode('utf-8'))

    # Get the relevant TOA object
    observation = get_object_or_404(models.Observation, pk=data['pk'])

    # Second, they have to belong to a group that has been granted edit privileges
    if not observation.can_edit(request.user):
        return HttpResponse(status=403)

    # Set the field value
    value = data['value']
    unit = u.Unit(data['unit'])
    if data['field'] in ["mjd_err"]:
        value = (float(value)*unit).to(observation_duration_units).value
    if data['field'] in ["freq", "bw"]:
        value = (float(value)*unit).to(observation_freq_units).value
    setattr(observation, data['field'], value)

    # Save the result to the database
    try:
        observation.save()
    except ValidationError as err:
        return HttpResponse(str(err), status=400)

    # Return "all is well"
    return HttpResponse(status=200)


@login_required
def lightcurve_view(request, pk):

    # Get the relevant Lightcurve object
    lightcurve = get_object_or_404(models.Lightcurve, pk=pk)

    # The user must be allowed to view this object
    if not lightcurve.can_view(request.user):
        return HttpResponse(status=403)

    # Get the (unique) working ephemeris (if any), so that any associated ToAs
    # can be worked out
    working_ephemeris = models.WorkingEphemeris.objects.filter(
        owner=request.user,
        ulp=lightcurve.ulp,
    ).first()

    context = {
        'lightcurve': lightcurve,
        'working_ephemeris': working_ephemeris,
    }

    lightcurve_points = lightcurve.points.all()

    if lightcurve_points.exists():

        # An attempt to make the plot client-side...
        pols = list({p.pol for p in lightcurve.points.all()})
        data = [
            {
                "x": list(lightcurve.times(pol=pol, dm=0.0)), # Setting dm=0 gets back un-dedispersed times
                "y": list(lightcurve.values(pol=pol)),
                "name": pol,
            } for pol in pols
        ]
        context['data'] = data

    return render(request, 'data/lightcurve.html', context)


@login_required
def lightcurve_add(request, pk):

    # Get the relevant Ulp object
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    if request.method == 'GET':

        context = {
            'ulp': ulp,
        }

        return render(request, 'data/lightcurve_new.html', context)

    elif request.method == 'POST':

        # First, try to open the file with np.loadtxt
        try:
            values = np.loadtxt(request.FILES['datafile'], ndmin=2)
        except:
            return HttpResponse("Could not open uploaded file with NumPy's loadtxt()", status=400)

        # Probably should clean/validate the data here...
        try:
            dm_freq = float(request.POST['dm_freq'])
        except:
            dm_freq = None

        # Add the lightcurve
        lightcurve = models.Lightcurve(
            owner=request.user,
            ulp=ulp,
            freq=request.POST['freq'],
            bw=request.POST['bw'],
            t0=request.POST['t0'],
            dt=request.POST['dt'],
            dm=request.POST['dm'],
            dm_freq=dm_freq,
            telescope=request.POST['telescope'],
        )
        lightcurve.save()

        # Add the lightcurve points
        pols = request.POST['pol_cols'].split()
        pols_without_underscore = [pol for pol in pols if pol != '_']
        no_repeats = (len(pols_without_underscore) == len(set(pols_without_underscore)))
        if not no_repeats:
            return HttpResponse(f"Cannot have repeated polarisations \"{request.POST['pol_cols']}\"", status=400)

        for p in range(len(pols)):
            pol = pols[p]

            # Ignore polarisations marked as '_'
            if pol == '_':
                continue

            # Add samples as LightcurvePoint objects
            for i in range(values.shape[0]):
                lightcurve_point = models.LightcurvePoint(
                    lightcurve=lightcurve,
                    sample_number=i,
                    pol=pol,
                    value=values[i,p],
                )
                # ^^^ YET TO DO: Figure out how to add errorbars!

                lightcurve_point.save()

        context = {
            'ulp': ulp,
        }

        return redirect('lightcurve_view', pk=lightcurve.pk)


@login_required
def pulsestack_view(request, pk):

    # Get the relevant Ulp object
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Get working ephemeris
    working_ephemeris = ulp.working_ephemerides.filter(owner=request.user).first()

    # If one doesn't exist for this user, make one!
    if not working_ephemeris:
        working_ephemeris = models.WorkingEphemeris(
            owner=request.user,
            ulp=ulp,
        )
        working_ephemeris.save()

    # Get all lightcurves that this owner can view
    lightcurves = permitted_to_view_filter(ulp.lightcurves.all(), request.user)

    # Extract the frequencies in order to build a colorscale
    freqs = [lc.freq for lc in lightcurves]
    freq_range = {'min': np.min(freqs), 'max': np.max(freqs)}

    # Get all the points, organised by lightcurve
    data = []
    for i in range(lightcurves.count()):
        lc = lightcurves[i]
        values = lc.values()
        values /= np.max(values)
        datum = {
            'idx': i,
            't0': lc.t0,
            'mjds': list(lc.bary_times(dm=0)),
            'values': list(values),
            'link': reverse('lightcurve_view', args=[lc.pk]),
            'freq_MHz': lc.freq,
        }
        data.append(datum)

    # Throw it all together into a context
    context = {
        'ulp': ulp,
        'working_ephemeris': working_ephemeris,
        'data': data,
        'freq_range': freq_range,
    }

    return render(request, 'data/pulsestack.html', context)


@login_required
def folding_view(request, pk):

    # Get the relevant Ulp object
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Get working ephemeris
    working_ephemeris = ulp.working_ephemerides.filter(owner=request.user).first()

    # If one doesn't exist for this user, make one!
    if not working_ephemeris:
        working_ephemeris = models.WorkingEphemeris(
            owner=request.user,
            ulp=ulp,
        )
        working_ephemeris.save()

    # Get all lightcurves that this owner can view
    lightcurves = permitted_to_view_filter(ulp.lightcurves.all(), request.user)

    # Get all the pulses from these light curves. Something will only show up on this
    # folding plot if there is an associated pulse object
    pulses = models.Pulse.objects.filter(lightcurve__in=lightcurves)

    # Extract the frequencies in order to build a colorscale
    freqs = [p.lightcurve.freq for p in pulses]
    freq_range = {'min': np.min(freqs), 'max': np.max(freqs)}

    data = []
    for i in range(pulses.count()):
        pulse = pulses[i]
        lc    = pulse.lightcurve
        pulse_times = barycentre(ulp, [pulse.mjd_start, pulse.mjd_end], sites[lc.telescope.lower()])
        mjd_ctr = (pulse_times[1] + pulse_times[0])/2
        date  = Time(mjd_ctr, scale='utc', format='mjd').isot

        # The peak_value_MJD, however, is not guaranteed to exist,
        # so must be dealt with separately to catch errors
        try:
            peak_value_MJD = barycentre(ulp, p.peak_value_MJD, sites[lc.telescope.lower()])
        except:
            peak_value_MJD = None

        datum = {
            'lc_pk': lc.pk,
            't0': lc.t0,
            'date': date,
            'telescope': lc.telescope,
            'link': reverse('lightcurve_view', args=[lc.pk]),
            'freq_MHz': lc.freq,
            'mjd_start': pulse_times[0],
            'mjd_ctr': mjd_ctr,
            'mjd_radius': (pulse_times[1] - pulse_times[0])/2,
            'mjd_end': pulse_times[1],
            'peak_value_MJD': peak_value_MJD or '',
            'peak_value_Jy': pulse.peak_value_Jy or '',
        }
        data.append(datum)

    # Throw it all together into a context
    context = {
        'ulp': ulp,
        'working_ephemeris': working_ephemeris,
        'data': data,
        'freq_range': freq_range,
    }

    return render(request, 'data/folding.html', context)


@login_required
def folding_toa_view(request, pk):

    # Get the relevant Ulp object
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Get working ephemeris
    working_ephemeris = ulp.working_ephemerides.filter(owner=request.user).first()

    # Get all ToAs for this ULP
    toas = models.Toa.objects.filter(template__working_ephemeris=working_ephemeris)

    def pack_data(toas, include_pk=True):
        return [
            {
                'toa_pk': toa.pk if include_pk else 0,
                'date': Time(toa.toa_mjd, scale='utc', format='mjd').isot,
                'link': reverse('toa_view', args=[toa.pk]) if include_pk else '',
                'toa_mjd': float(toa.toa_mjd),
                'toa_err_s': toa.toa_err_s,
                'ampl': toa.ampl,
                'ampl_err': toa.ampl_err,
                'freq_MHz': toa.pulse.lightcurve.freq,
                'include_in_fit': 'true' if toa.include_in_fit else 'false',
            } for toa in toas
        ]
    data = pack_data(toas)

    # Throw it all together into a context
    context = {
        'toas': toas,
        'ulp': ulp,
        'working_ephemeris': working_ephemeris,
        'data': data,
        #'fit': {
        #   'working_ephemeris': fitted_working_ephemeris,
        #    'data': predicted_data,
        #},
    }

    return render(request, 'data/folding_toa.html', context)


@login_required
def update_working_ephemeris(request, pk):

    # Get the relevant WorkingEphemeris object
    working_ephemeris = get_object_or_404(models.WorkingEphemeris, pk=pk)

    # If this is a POST request, save the provided values, also making sure the user is allowed to edit it
    if request.method == "POST" and working_ephemeris.can_edit(request.user):

        # Do some validation here...

        for field in ['ra', 'dec', 'pepoch', 'p0', 'p1', 'pb', 'dm', 'tausc_1GHz', 'spec_alpha', 'spec_q']:
            try:
                setattr(working_ephemeris, field, float(request.POST[field]))
            except:
                pass

        # Now update the covariance, if there is anything to update
        # First, drop the old one, if it exists
        if request.POST.get('covariance'):
            json_str = ["'" + request.POST['covariance'] + "'"]
            covariance = serializers.WorkingEphemerisCovarianceDeserializer(json_str)[0].object
            covariance.save()
            if not working_ephemeris.covariance:
                working_ephemeris.covariance = covariance

        working_ephemeris.save()

    next = request.POST.get('next', reverse('timing_choose_ulp'))

    return HttpResponseRedirect(next)


@login_required
def update_selected_working_ephemeris(request):

    # Turn the data into a dictionary
    data = json.loads(request.body.decode('utf-8'))

    # Get the relevant WorkingEphemeris object
    working_ephemeris = models.WorkingEphemeris.objects.get(pk=data['pk'])
    if working_ephemeris is None:
        return JsonResponse({'message': f"No working ephemeris found with pk = {data['pk']}"}, status=400)

    # If this is a POST request, save the provided values, also making sure the user is allowed to edit it
    if request.method == "PUT" and working_ephemeris.can_edit(request.user):

        new_ephemeris_values = {} # This is for returning to the client so that webpage values can be updated

        for field in ['ra', 'dec', 'pepoch', 'p0', 'p1', 'dm', 'p_aw', 't0_aw', 'duration_aw']:
            try:
                if data[f'select_{field}'] == False:
                    continue
            except:
                continue

            try:
                value = float(data[field])
            except:
                value = data[field]

            new_ephemeris_values[field] = value
            # This ^^^ makes sure the page gets populated with the original value,
            # even if the corresponding model isn't updated

            try:
                setattr(working_ephemeris, field, value)
            except:
                pass

        working_ephemeris.save()

        new_ephemeris_values['pk'] = data['pk']

        # Now update the (whole) covariance matrix.
        # First, deserialize it
        covariance = list(serializers.WorkingEphemerisCovarianceDeserializer(data['covariance']))[0].object

        # Ensure that this covariance is being applied to the correct working ephemeris
        # by forcing
        if working_ephemeris.covariance:
            working_ephemeris.covariance.delete()

        working_ephemeris.covariance = covariance
        covariance.save()
        working_ephemeris.save()

        return JsonResponse(new_ephemeris_values, status=200)

    return JsonResponse({'message': f"User does not have edit privileges for this ephemeris"}, status=400)


@login_required
def toa_view(request, pk):

    # Get the relevant ToA object
    toa = get_object_or_404(models.Toa, pk=pk)

    sc_idx = request.GET.get("sc_idx")
    try:
        sc_idx = float(sc_idx)
    except:
        sc_idx = -4.0

    # Get lightcurve shorthand
    lc = toa.pulse.lightcurve
    times = lc.bary_times(dm=0.0)
    values = lc.values()

    min_t, max_t = toa.get_toa_range()

    N = 500
    template_times = np.linspace(min_t, max_t, N)
    template_values = toa.ampl * toa.template.values(template_times - float(toa.toa_mjd))

    # If there is a scattering timescale associated with this ToA, show the scattered
    # template as well
    try:
        tausc_1GHz = float(request.GET.get('tausc_1GHz'))
    except:
        tausc_1GHz = toa.template.working_ephemeris.tausc_1GHz

    if tausc_1GHz is not None:
        tausc = tausc_1GHz * (lc.freq/1e3)**sc_idx
        scattered_template_values = toa.ampl * toa.template.values(
                template_times - float(toa.toa_mjd),
                tausc=tausc,
                freq=lc.freq,
                bw=lc.bw,
                sc_idx=sc_idx)
    else:
        scattered_template_values = None

    # Deal with the baseline
    baseline_degree = -1 # Default means "no baseline fitting was used"

    if toa.baseline_level is not None:
        template_values += toa.baseline_level
        if scattered_template_values is not None:
            scattered_template_values += toa.baseline_level
        baseline_degree = 0 # Constant function

    if toa.baseline_slope is not None:
        slope_adjustment = toa.baseline_slope*(template_times - float(toa.toa_mjd))
        template_values += slope_adjustment
        if scattered_template_values is not None:
            scattered_template_values += slope_adjustment
        baseline_degree = 1 # Linear function

    toa_err_days = toa.toa_err_s / 86400.0

    context = {
        'toa': toa,
        'toa_err_days': toa_err_days,
        'times': list(times),
        'values': list(values),
        'template_times': list(template_times),
        'template_values': list(template_values),
        'scattered_template_values': list(scattered_template_values) if scattered_template_values is not None else None,
        'tausc_1GHz': tausc_1GHz or toa.template.working_ephemeris.tausc_1GHz,
        'baseline_degree': baseline_degree,
    }

    return render(request, 'data/toa.html', context)


@login_required
def refit_toa(request, pk):

    # Get the relevant Toa object
    toa = get_object_or_404(models.Toa, pk=pk)

    sc_idx = request.GET.get("sc_idx")
    try:
        sc_idx = float(sc_idx)
    except:
        sc_idx = -4.0

    # Get fitting options
    baseline_degree = int(request.POST.get('baseline_degree')) # Only values of 0 and 1 currently carry meaning. Everything else means "don't fit"
    use_tausc_1GHz = request.POST.get('use_tausc_1GHz') == 'on'
    try:
        tausc_1GHz = float(request.POST.get('tausc_1GHz'))
        tausc = tausc_1GHz * (toa.pulse.lightcurve.freq/1e3)**sc_idx if use_tausc_1GHz else None
    except:
        tausc = None

    # Refit the ToA
    if baseline_degree == -1:
        toa.refit(baseline_level=None, baseline_slope=None, tausc=tausc)
    elif baseline_degree == 0:
        toa.refit(baseline_level=0.0, baseline_slope=None, tausc=tausc)
    elif baseline_degree == 1:
        toa.refit(baseline_level=0.0, baseline_slope=0.0, tausc=tausc)
    else:
        toa.refit(tausc=tausc)

    url = reverse('toa_view', args=[pk])
    query_string = urlencode({'tausc_1GHz': tausc_1GHz}) if use_tausc_1GHz else toa.template.working_ephemeris.tausc_1GHz
    url_with_query = f"{url}?{query_string}"

    return redirect(url_with_query)


@login_required
def toa_for_pulse(request, pk):

    # Get the relevant Lightcurve object
    pulse = get_object_or_404(models.Pulse, pk=pk)

    # There must be an existing (unique) working ephemeris
    we = get_object_or_404(models.WorkingEphemeris, owner=request.user, ulp=pulse.lightcurve.ulp)

    # There must also be an existing template
    # At the moment, this will fail if multiple templates exist that fit these criteria,
    # which could happen because (AFAIK) there isn't any such constraint on the database table
    # (but perhaps there should be...??)
    template = models.Template.objects.filter(owner=request.user, working_ephemeris=we).first()

    # If a ToA already exists with this pulse number, just go to that page;
    # otherwise, create one
    toa = models.Toa.objects.filter(pulse=pulse, template=template).first()

    if toa is None:
        print("here3.0")
        # Make an initial fit with the default settings
        toa = models.Toa(
            pulse=pulse,
            template=template,
        )

        toa.refit(toa_mjd='peak', ampl='peak')
        # ^^^ This saves the new Toa object by default

    # Navigate to the ToA page
    return redirect('toa_view', pk=toa.pk)


@login_required
def download_toas(request, pk):

    # Get the relevant Ulp object
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Get working ephemeris
    working_ephemeris = ulp.working_ephemerides.filter(owner=request.user).first()

    # Get all ToAs for this ULP
    toas = models.Toa.objects.filter(template__working_ephemeris=working_ephemeris)

    current_datetime = now()
    current_datetime_pretty = current_datetime.strftime('%Y-%m-%d %H:%M:%S')
    current_datetime_filename = current_datetime.strftime('%Y%m%d')
    ulp_name_filename = ulp.name.replace(' ', '_')
    ascii_content = "MODE 1\nFORMAT 1"
    ascii_content += f"\nC Generated by the ULP database {current_datetime}"
    ascii_content += f"\nC Timing residuals for {ulp.name}"
    for toa in toas:
        lc = toa.pulse.lightcurve
        telescope = lc.telescope
        freq = lc.freq
        # vvv Elaborate shenanigans to get BACK to topocentric ToAs
        mjd = toa.toa_mjd + Decimal(lc.t0 - barycentre(ulp, [lc.t0], location=sites[lc.telescope.lower()])[0])
        mjd_err_us = toa.toa_err_s * 1e6
        backend = telescope
        ascii_content += f"\n{toa.pk} {freq:f} {mjd} {mjd_err_us} {backend}"

    response = HttpResponse(ascii_content, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename="{ulp_name_filename}_{current_datetime_filename}.tim"'

    # Return the response!
    return response


def download_times_of_arrival(request, pk):

    # Retrieve the selected ToAs from the specified ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)
    toas = ulp.times_of_arrival.all()
    # Now limit only to those that this user has view access to
    toas = permitted_to_view_filter(toas, request.user)

    only_published = request.GET.get('only_published') == '1'
    if only_published:
        toas = toas.filter(published_in__isnull=False)

    current_datetime = now()
    current_datetime_pretty = current_datetime.strftime('%Y-%m-%d %H:%M:%S')
    current_datetime_filename = current_datetime.strftime('%Y%m%d')
    ulp_name_filename = ulp.name.replace(' ', '_')
    ascii_content = "MODE 1\nFORMAT 1"
    ascii_content += f"\nC Generated by the ULP database {current_datetime}"
    ascii_content += f"\nC Timing residuals for {ulp.name}"
    for toa in toas:
        if toa.published_in is not None:
            name = re.sub(r'[^a-zA-Z0-9.]', '_', toa.published_in.citet_text)
            name = re.sub(r'_+', '_', name)
            name = name.strip('_')
        else:
            name = f'{toa.pk}'
        mjd = toa.raw_mjd
        freq = toa.freq
        mjd_err_us = float(toa.mjd_err) * 86400e6
        backend = toa.telescope_name.replace(" ", "_")
        ascii_content += f"\n{name} {freq:f} {mjd} {mjd_err_us} {backend}"

    response = HttpResponse(ascii_content, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename="{ulp_name_filename}_{current_datetime_filename}.tim"'

    # Return the response!
    return response


@login_required
def download_working_ephemeris(request, pk):

    # Get the relevant WorkingEphemeris object (which the user must own)
    working_ephemeris = get_object_or_404(models.WorkingEphemeris, pk=pk, owner=request.user)

    # Get/generate some values to be placed in the ephemeris
    ulp = working_ephemeris.ulp
    ulp_name_filename = ulp.name.replace(' ', '_')
    ascii_content = f"PSRJ            {ulp_name_filename}"

    ra = ulp.measurements.filter(parameter__name="Right ascension").first()
    if ra is not None:
        ascii_content += f"\nRAJ             {ra.formatted_quantity}"

    dec = ulp.measurements.filter(parameter__name="Declination").first()
    if dec is not None:
        ascii_content += f"\nDECJ            {dec.formatted_quantity}"

    ascii_content += f"\nDM              {working_ephemeris.dm}"
    ascii_content += f"\nPEPOCH          {working_ephemeris.pepoch}"
    ascii_content += f"\nPOSEPOCH        {working_ephemeris.pepoch}"
    ascii_content += f"\nDMEPOCH         {working_ephemeris.pepoch}"
    ascii_content += f"\nF0              {1/working_ephemeris.p0}"
    ascii_content += f"\nF1              {-working_ephemeris.p1/working_ephemeris.p0**2}"
    ascii_content += f"\nEPHVER          TEMPO2"
    ascii_content += f"\nUNITS           TCB"

    current_datetime_filename = now().strftime('%Y%m%d')
    response = HttpResponse(ascii_content, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename="{ulp_name_filename}_{current_datetime_filename}.par"'

    # Return the response!
    return response


@api_view(['POST'])
@authentication_classes([TokenAuthentication])
@permission_classes([IsAuthenticated])
def upload_lightcurve(request):
    if 'file' in request.FILES:
        uploaded_file = request.FILES['file']

        try:
            file_content = uploaded_file.read()
            stream = io.BytesIO(file_content)
            data = np.load(stream, allow_pickle=True).item() # This should now be a dictionary

            ###########
            try:
                ulp_pk = request.GET.get('ulp')
                if not ulp_pk or not ulp_pk.isdigit():
                    raise ValueError(f"Invalid ulp id {ulp_pk}")
                ulp = published_models.Ulp.objects.get(pk=int(ulp_pk))
            except published_models.Ulp.DoesNotExist:
                raise Exception(f"ULP {ulp_pk} does not exist.")
            except ValueError as e:
                raise Exception(f"Invalid ulp id: {ulp_pk}. Error: {str(e)}")
            except Exception as e:
                raise Exception(f"An unexpected error occurred: {str(e)}")
            
            # Add the lightcurve
            lightcurve = models.Lightcurve(
                owner=request.user,
                ulp=ulp,
                freq=data['CTR_FREQ'],
                bw=data['BW'],
                t0=data['TIME'][0],
                dt=(data['TIME'][-1] - data['TIME'][0])/(len(data['TIME']) - 1) * 86400.0,
                dm=data['DM'],
                dm_freq=data['DMREFFREQ'],
                telescope=data['TELESCOPE'],
            )

            # If the DM frequency is infinite, store 'None' in the database, which
            # has that exact special meaning
            if lightcurve.dm_freq == np.inf:
                lightcurve.dm_freq = None

            lightcurve.save()

            # Add the lightcurve points
            for i in range(data['LIGHTCURVE'].shape[0]):
                lightcurve_point = models.LightcurvePoint(
                    lightcurve=lightcurve,
                    sample_number=i,
                    pol=data['POL'],
                    value=data['LIGHTCURVE'][i],
                )

                lightcurve_point.save()

            ###########

            return JsonResponse({'message': 'Lightcurve uploaded successfully', 'pk': lightcurve.pk, 'url': reverse('lightcurve_view', args=[lightcurve.pk])}, status=201)
        except Exception as e:
            return JsonResponse({'error': 'Failed to process the pickle file', 'details': str(e)}, status=400)

    return JsonResponse({'error': 'Invalid request'}, status=400)


@authentication_classes([TokenAuthentication])
@permission_classes([IsAuthenticated])
def fit_ephemeris(request, ulp_pk):

    data = json.loads(request.body.decode('utf-8'))

    # Retrieve the selected ULP and working ephemeris
    ulp = get_object_or_404(published_models.Ulp, pk=ulp_pk)

    we = models.WorkingEphemeris(
        owner=request.user,
        ulp=ulp,
        ra=data['ra'],
        dec=data['dec'],
        pepoch=float(data['pepoch']),
        p0=float(data['p0']),
        p1=float(data['p1']),
        dm=float(data['dm']),
    )

    toas_data = toa_data(request.user, we)

    bounds_dict = {
        'ra': [0.0, 360.0],
        'dec': [-90.0, 90.0],
        'pepoch': [-np.inf, np.inf],
        'p0': [0, np.inf],
        'p1': [-np.inf, np.inf],
        'dm': [0, np.inf],
    }

    # Initialise the answer to be the same thing as the input
    best_fit_ephemeris = {}
    for param in ['ra', 'dec', 'pepoch', 'p0', 'p1', 'dm']:
        best_fit_ephemeris[param] = data[param]

    # RA-No  DEC-No  PEPOCH-Yes  P0-Yes  P1-No  DM-No
    if not data['select_ra'] and not data['select_dec'] and data['select_pepoch'] and data['select_p0'] and not data['select_p1'] and not data['select_dm']:
        init_pepoch = Time(data['pepoch'], scale='utc', format='mjd')
        def func(toas, pepoch, p0_s):
            # In order to avoid convergence issues, toas and pepoch in this function are relative
            # to the initial guess pepoch
            toas_mjd = toas + init_pepoch.mjd
            pepoch_mjd = pepoch + init_pepoch.mjd
            ephemeris = {'pepoch': pepoch_mjd, 'p0': p0_s, 'p1': data['p1']}
            time = Time(toas_mjd, scale='utc', format='mjd')
            pulse_phases = calc_pulse_phase(time, ephemeris)
            nearest_pulse_phases = np.round(pulse_phases)
            predicted_mjds = calc_mjd(nearest_pulse_phases, ephemeris).mjd
            return predicted_mjds - init_pepoch.mjd

        x = [toa['mjd'] - calc_dmdelay(we.dm*u.pc/u.cm**3, toa['freq_MHz']*u.MHz, np.inf*u.MHz).to('d').value + toa['bc_correction'] - init_pepoch.mjd for toa in toas_data]
        sigma = [toa['mjd_err'] for toa in toas_data]
        p0 = [0.0, data['p0']]
        bounds = [
            (bounds_dict['pepoch'][0], bounds_dict['p0'][0]),
            (bounds_dict['pepoch'][1], bounds_dict['p0'][1]),
        ]

        popt, pcov = curve_fit(func, x, x, p0=p0, bounds=bounds, sigma=sigma)

        best_fit_ephemeris['pepoch'] = f"{popt[0] + init_pepoch.mjd}"
        best_fit_ephemeris['p0'] = f"{popt[1]}"

        # Populate a covariance
        covariance = models.WorkingEphemerisCovariance(
            pepoch_pepoch=pcov[0,0],
            pepoch_p0=pcov[0,1],
            p0_p0=pcov[1,1],
        )

    elif not data['select_ra'] and not data['select_dec'] and data['select_pepoch'] and data['select_p0'] and data['select_p1'] and not data['select_dm']:
        # RA-No  DEC-No  PEPOCH-Yes  P0-Yes  P1-Yes  DM-No
        p1_scale = 1e12
        def func(toas_mjd, pepoch_mjd, p0_s, p1__scaled):
            ephemeris = {'pepoch': pepoch_mjd, 'p0': p0_s, 'p1': p1__scaled / p1_scale}
            time = Time(toas_mjd, scale='utc', format='mjd')
            pulse_phases = calc_pulse_phase(time, ephemeris)
            nearest_pulse_phases = np.round(pulse_phases)
            predicted_mjds = calc_mjd(nearest_pulse_phases, ephemeris).mjd
            return predicted_mjds

        x = [toa['mjd'] - calc_dmdelay(we.dm*u.pc/u.cm**3, toa['freq_MHz']*u.MHz, np.inf*u.MHz).to('d').value + toa['bc_correction'] for toa in toas_data]
        sigma = [toa['mjd_err'] for toa in toas_data]
        p0 = [data['pepoch'], data['p0'], u.Quantity(data['p1'])*p1_scale]
        bounds = [
            (bounds_dict['pepoch'][0], bounds_dict['p0'][0], bounds_dict['p1'][0]),
            (bounds_dict['pepoch'][1], bounds_dict['p0'][1], bounds_dict['p1'][1]),
        ]

        popt, pcov = curve_fit(func, x, x, p0=p0, bounds=bounds, sigma=sigma)

        best_fit_ephemeris['pepoch'] = f"{popt[0]}"
        best_fit_ephemeris['p0'] = f"{popt[1]}"
        best_fit_ephemeris['p1'] = f"{popt[2]/p1_scale}"

        # Populate a covariance
        covariance = models.WorkingEphemerisCovariance(
            pepoch_pepoch=pcov[0,0],
            pepoch_p0=pcov[0,1],
            pepoch_p1=pcov[0,2]/p1_scale,
            p0_p0=pcov[1,1],
            p0_p1=pcov[1,2]/p1_scale,
            p1_p1=pcov[2,2]/p1_scale**2,
        )

    else:
        return JsonResponse({'message': 'The chosen combination of fit parameters is not yet supported'}, status=400)

    # Serialize the covariance matrix
    cov = serializers.WorkingEphemerisCovarianceSerializer().serialize([covariance])

    return JsonResponse([best_fit_ephemeris, cov], safe=False, status=200)



@api_view(['POST'])
@authentication_classes([TokenAuthentication])
@permission_classes([IsAuthenticated])
def write_toas(request):

    if 'toa_file' not in request.FILES:
        raise rest_exceptions.ParseError("No file uploaded")

    toa_file = request.FILES['toa_file']

    fmt = request.data.get('format', 'tim_format_1')

    # Get overwrite/add/ignore mode
    mode = request.data.get('mode')
    supported_modes = ['overwrite', 'add', 'ignore']
    if mode not in supported_modes:
        raise rest_exceptions.ParseError(f"Mode {mode} not supported. Must be one of {supported_modes}.")

    # Get specified LPT (required)
    ulp = published_models.Ulp.objects.filter(name=request.data.get('lpt')).first()
    if ulp is None:
        raise rest_exceptions.ParseError(f"Unrecognised LPT: '{request.data.get('lpt')}'")

    # Get user's working ephemeris (required)
    we = models.WorkingEphemeris.objects.filter(owner=request.user, ulp=ulp).first()
    if we is None:
        raise rest_exceptions.ParseError(f"User {request.user} has no working ephemeris for {ulp}")

    supported_formats = ['astropy_qtable', 'tim_format_1']
    if fmt not in supported_formats:
        raise rest_exceptions.ParseError(f"{fmt} is not a supported format")

    if fmt.lower() == 'astropy_qtable':
        try:
            table = QTable.read(toa_file, format='ascii.ecsv')
        except:
            raise rest_exceptions.ParseError('The uploaded file could not be parsed as an AstroPy QTable')

        # Make sure the necessary columns are present
        columns = table.keys()
        required_columns = ['ToA', 'ToA_err', 'freq', 'telescope']
        for required_column in required_columns:
            if required_column not in columns:
                raise rest_exception.ParseError(f"QTable missing '{required_column}' column")

        # Generate list of ToAs based on table contents
        toas = [
            models.TimeOfArrival(
                owner=request.user,
                ulp=ulp,
                raw_mjd=row.get('ToA'),
                mjd=row.get('ToA'),
                mjd_err=row.get('ToA_err').to('d').value,
                telescope_name=row.get('telescope'),
                freq=row.get('freq').to('MHz').value,
                fluence_Jy_s=row.get('fluence').to('Jy s').value if row.get('fluence') is not None else None,
                peak_flux_Jy=row.get('fitted_peak_flux_density').to('Jy').value if row.get('fitted_peak_flux_density') is not None else None,
                barycentred=False,
                dedispersed=False,
            ) for row in table
        ]

    elif fmt.lower() == 'tim_format_1':
        lines = toa_file.read().decode().splitlines()

        # Assume the first two lines ("FORMAT 1\nMODE 1") are required by the format.
        # Also, ignore any line starting with "C", which indicates a comment.
        good_lines = [line.split() for line in lines[2:] if line[0] != "C" and len(line) > 0]

        # Build ToAs out of the good lines
        toas = [
            models.TimeOfArrival(
                owner=request.user,
                ulp=ulp,
                raw_mjd=float(line[2]),
                mjd=float(line[2]),
                mjd_err=(float(line[3])*u.us).to('d').value,
                telescope_name=line[4],
                freq=float(line[1]),
            ) for line in good_lines
        ]

    else:
        raise rest_exceptions.ParseError(f"Format '{fmt}' not yet implemented. Working on it!")

    for toa in toas:
        matching_toa = permitted_to_edit_filter(models.TimeOfArrival.objects.filter(
            ulp=ulp,
            freq__gte=toa.freq*0.98,
            freq__lte=toa.freq*1.02,
            raw_mjd__gte=toa.raw_mjd - we.p0/2/86400,
            raw_mjd__lte=toa.raw_mjd + we.p0/2/86400,
            telescope_name=toa.telescope_name,
        ), request.user).first()

        if matching_toa is not None: # i.e. there is an existing ToA with same Telescope at same frequency catching the same pulse
            print(f"Found matching ToA: {matching_toa}")
            if mode == 'ignore':
                continue
            elif mode == 'add':
                toa.save()
            elif mode == 'overwrite':
                matching_toa.raw_mjd = toa.raw_mjd
                matching_toa.mjd = toa.mjd
                matching_toa.mjd_err = toa.mjd_err
                matching_toa.freq = toa.freq
                matching_toa.fluence_Jy_s = toa.fluence_Jy_s
                matching_toa.peak_flux_Jy = toa.peak_flux_Jy
                matching_toa.barycentred = toa.barycentred
                matching_toa.dedispersed = toa.dedispersed
                matching_toa.save()
            else:
                # Should never get here: This check is done earlier in this function.
                raise rest_exceptions.ParseError(f"Mode {mode} not supported. Must be one of {supported_modes}.")
        else:
            toa.save()

    return JsonResponse('Success', safe=False, status=200)
