from django.shortcuts import render, get_object_or_404, redirect
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.db.models import Q, Value, BooleanField
from django.urls import reverse
from django.core.exceptions import ValidationError
from . import models
from common.utils import *
from django.utils.timezone import now
from urllib.parse import urlencode

from rest_framework.authentication import TokenAuthentication
from rest_framework.permissions import IsAuthenticated
from rest_framework.decorators import api_view, authentication_classes, permission_classes

from published import models as published_models
from published.views import get_accessible_measurements

import numpy as np
from scipy.optimize import curve_fit
import json
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz, get_sun
from astropy.constants import c
from astropy.time import Time
from decimal import Decimal
import io

from spiceypy.spiceypy import spkezr, furnsh, j2000, spd, unload
import pandas as pd
from plotly.offline import plot
import plotly.express as px
#import plotly.graph_objects as go
#from plotly.subplots import make_subplots

site_names = sorted(EarthLocation.get_site_names(), key=str.casefold)

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

def bc_corr(coord, times, ephemeris_file='de430.bsp'):
    '''
    coord = SkyCoord object (from astropy) representing the location of the source
    times = Time array of MJDs
    ephemeris_file = e.g. de430.bsp
    '''
    try:
        furnsh(ephemeris_file)
    except:
        raise Exception("Cannot load ephemeris file {}\n".format(ephemeris_file))
    #jds = times.mjd + 2400000.5
    ets = (times.jd - j2000())*spd()
    r_earths = [spkezr("earth", et, "j2000", "NONE", "solar system barycenter")[0][:3] for et in ets]
    r_src_normalised = [
        np.cos(coord.ra.rad)*np.cos(coord.dec.rad),
        np.sin(coord.ra.rad)*np.cos(coord.dec.rad),
        np.sin(coord.dec.rad),
    ]
    delays = np.array([np.dot(r_earth, r_src_normalised) for r_earth in r_earths]) * u.km / c # (spkezr returns km)

    return delays


def calc_pulse_phase(time, ephemeris):
    '''
    time is an astropy Time object
    ephemeris['pepoch'] should be in days
    ephemeris['p0'] should be in seconds
    '''
    try:
        pepoch = Time(ephemeris['pepoch'], format='mjd')
        p0 = ephemeris['p0']*u.s
        return ((time - pepoch)/p0).decompose()
    except:
        return 0.0


def calc_mjd(pulse_phase, ephemeris):
    '''
    This is the inverse of calc_pulse_phase()
    '''
    pepoch = Time(ephemeris['pepoch'], format='mjd')
    p0 = ephemeris['p0']*u.s
    return pepoch + pulse_phase*p0


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
    mjds = Time([float(toa.raw_mjd) for toa in toas], format='mjd')
    bc_corrections = bc_corr(we.coord, mjds).to('day').value

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

    # Barycentre
    mjds = Time([float(obs.start_mjd) for obs in obss], format='mjd') + ([obs.duration for obs in obss] * u.s)/2
    bc_corrections = bc_corr(we.coord, mjds).to('day').value

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


@login_required
def timing_choose_ulp_view(request):

    # Get ULPs to which they have access
    ulps = published_models.Ulp.objects.filter(
        Q(data_access_groups__user=request.user) |
        Q(whitelist_users=request.user)
    ).distinct()

    context = {
        'ulps': ulps,
    }

    return render(request, 'data/timing_choose_ulp.html', context)


def get_toa_predictions(start, end, freq_MHz, pepoch, p0, dm, telescope, coord, output_toa_format='mjd'):

    # First, assume the given mjd_start and mjd_end are in fact topocentric dispersed MJDs,
    # so to get the right range, convert them to dedispersed, barycentric
    time_range = Time([start, end])
    dmdelay = calc_dmdelay(dm, freq_MHz, np.inf*u.MHz)
    time_range -= dmdelay
    time_range += bc_corr(coord, time_range)

    predicted_barycentric_toas = generate_toas(time_range[0], time_range[1], ephemeris)
    predicted_topocentric_toas = predicted_barycentric_toas - bc_corr(coord, predicted_barycentric_toas)
    predicted_dispersed_toas = predicted_topocentric_toas + dmdelay

    altaz = coord.transform_to(AltAz(obstime=predicted_dispersed_toas, location=EarthLocation.of_site(telescope)))
    sun_altaz = [get_sun(predicted_dispersed_toas[i]).transform_to(AltAz(obstime=predicted_dispersed_toas[i], location=EarthLocation.of_site(telescope))) for i in range(len(predicted_barycentric_toas))]

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
        } for i in range(len(predicted_barycentric_toas)) if altaz[i].alt.value >= min_el and sun_altaz[i].alt.value < max_sun_el
    ]

    return predicted_toas


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

    # Get the working ephemerides that are available to this user
    working_ephemerides = permitted_to_view_filter(models.WorkingEphemeris.objects.filter(ulp=ulp), request.user).order_by('owner')

    # If there aren't any, make a default one for the user!
    if not working_ephemerides.filter(owner=request.user).exists():
        we = models.WorkingEphemeris(owner=request.user, ulp=ulp)
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
    context['selected_working_ephemeris'] = selected_working_ephemeris

    # Pool together lists of published values to offer as options to the user
    context['published_ephemeris_values'] = {
        'ra': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Right ascension", article__isnull=False),
        'dec': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Declination", article__isnull=False),
        'pepoch': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="PEPOCH", article__isnull=False),
        'p0': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Period", article__isnull=False),
        'dm': published_models.Measurement.objects.filter(ulp=ulp, parameter__name="Dispersion measure", article__isnull=False),
    }

    # Construct a dictionary out of the ephemeris
    ephemeris = {
        'ra': selected_working_ephemeris.ra,
        'dec': selected_working_ephemeris.dec,
        'pepoch': selected_working_ephemeris.pepoch,
        'p0': selected_working_ephemeris.p0,
        'dm': selected_working_ephemeris.dm,
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
        dm = float(request.POST.get('dm', '0.0'))
        ra = request.POST.get('ra', '00:00:00')
        dec = request.POST.get('dec', '00:00:00')

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
        ephemeris['dm'] = dm
        ephemeris['ra'] = ra
        ephemeris['dec'] = dec

        # If they've also provided other form values, make a table of predicted values
        if mjd_start is not None and mjd_end is not None and mjd_dispersion_frequency is not None and pepoch is not None and p0 is not None and telescope is not None:
            context['predicted_toas'] = get_toa_predictions(
                mjd_range[0],                       # Start of time range
                mjd_range[-1],                      # End of time range
                mjd_dispersion_frequency*u.MHz,     # Observing frequency
                pepoch,                             # PEPOCH (MJD)
                p0,                                 # Period
                ephemeris['dm']*u.pc/u.cm**3,       # DM
                telescope,                          # Telescope string
                coord,                              # Target source coordinates
                output_toa_format=output_toa_format # Desired output ToA format
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
        pulse_times = barycentre(ulp, [pulse.mjd_start, pulse.mjd_end], EarthLocation.of_site(lc.telescope))
        mjd_ctr = (pulse_times[1] + pulse_times[0])/2
        date  = Time(mjd_ctr, scale='utc', format='mjd').isot

        # The peak_value_MJD, however, is not guaranteed to exist,
        # so must be dealt with separately to catch errors
        try:
            peak_value_MJD = barycentre(ulp, p.peak_value_MJD, EarthLocation.of_site(lc.telescope))
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

    # How to get residuals for an arbitrary working ephemeris + toas:
    #residuals_for_toa_qs(self, toa_qs)
    '''
    # Make a new working ephemeris to hold the best fitting solution
    def fold_for_curve_fit(pulses, pepoch, period):
        # Use the folding function defined in the WorkingEphemeris model
        # to produce the predicted phases
        we = models.WorkingEphemeris(
            ulp=ulp,
            pepoch=pepoch,
            p0=period,
        )
        mjds = we.unfold(pulse_numbers, np.zeros(pulse_numbers.shape))

        return mjds

    x = np.array([toa.pulse_number for toa in toas if toa.include_in_fit])
    y = np.array([toa.toa_mjd for toa in toas if toa.include_in_fit])
    p0 = [working_ephemeris.pepoch, working_ephemeris.p0] # 1st "p0" means "initial parameters for curve_fit"; 2nd "p0" means "period"
    popt, pcov = curve_fit(fold_for_curve_fit, x, y, p0=p0)
    covariance = models.WorkingEphemerisCovariance(
        pepoch_pepoch=pcov[0,0],
        pepoch_p0=pcov[0,1],
        p0_p0=pcov[1,1],
    )
    predicted_y = fold_for_curve_fit(x, *popt)
    fitted_working_ephemeris = models.WorkingEphemeris(
        ulp=ulp,
        pepoch=popt[0],
        p0=popt[1],
        covariance=covariance,
    )

    # Create predicted ToAs based on fit
    predicted_toas = [
        models.Toa(
            pulse_number=x[i],
            template=toas[i].template,
            toa_mjd=Decimal(predicted_y[i]),
            toa_err_s=0.0,
            ampl=0.0,
            ampl_err=0.0,
            ampl_ref_freq=0.0,
        ) for i in range(len(x))
    ]
    predicted_data = pack_data(predicted_toas, include_pk=False)
    '''

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
        covariance = models.WorkingEphemerisCovariance.from_str(request.POST['covariance_list'])

        if covariance is not None:
            if working_ephemeris.covariance is not None:
                working_ephemeris.covariance.delete()
            covariance.save()
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

        for field in ['ra', 'dec', 'pepoch', 'p0', 'dm']:
            try:
                value = float(data[field])
            except:
                value = data[field]

            new_ephemeris_values[field] = value
            # This ^^^ makes sure the page gets populated with the original value,
            # even if the corresponding model isn't updated

            try:
                if data[f'select_{field}'] == False:
                    continue
            except:
                pass

            try:
                setattr(working_ephemeris, field, value)
            except:
                pass

        working_ephemeris.save()

        new_ephemeris_values['pk'] = data['pk']
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
        mjd = toa.toa_mjd + Decimal(lc.t0 - barycentre(ulp, [lc.t0], location=EarthLocation.of_site(lc.telescope))[0])
        mjd_err_us = toa.toa_err_s * 1e6
        backend = telescope
        ascii_content += f"\n{toa.pk} {freq:f} {mjd} {mjd_err_us} {backend}"

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
    ascii_content += f"\nF1              {working_ephemeris.p1/working_ephemeris.p0**2}"
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
        dm=float(data['dm']),
    )

    toas_data = toa_data(request.user, we)

    bounds_dict = {
        'ra': [0.0, 360.0],
        'dec': [-90.0, 90.0],
        'pepoch': [-np.inf, np.inf],
        'p0': [0, np.inf],
        'dm': [0, np.inf],
    }

    # Initialise the answer to be the same thing as the input
    best_fit_ephemeris = {}
    for param in ['ra', 'dec', 'pepoch', 'p0', 'dm']:
        best_fit_ephemeris[param] = data[param]

    # RA-No  DEC-No  PEPOCH-Yes  P0-Yes  DM-No
    if not data['select_ra'] and not data['select_dec'] and data['select_pepoch'] and data['select_p0'] and not data['select_dm']:
        func = fit_ephemeris_pepoch_p0
        x = [toa['mjd'] - calc_dmdelay(we.dm*u.pc/u.cm**3, toa['freq_MHz']*u.MHz, np.inf*u.MHz).to('d').value + toa['bc_correction'] for toa in toas_data]
        sigma = [toa['mjd_err'] for toa in toas_data]
        p0 = [data['pepoch'], data['p0']]
        bounds = [
            (bounds_dict['pepoch'][0], bounds_dict['p0'][0]),
            (bounds_dict['pepoch'][1], bounds_dict['p0'][1]),
        ]

        popt, pcov = curve_fit(func, x, x, p0=p0, bounds=bounds, sigma=sigma)

        best_fit_ephemeris['pepoch'] = f"{popt[0]}"
        best_fit_ephemeris['p0'] = f"{popt[1]}"

        # Populate a covariance
        cov = {
            'pepoch_pepoch': pcov[0,0],
            'pepoch_p0': pcov[0,1],
            'p0_p0': pcov[1,1],
        }

    else:
        return JsonResponse({'message': 'The chosen combination of fit parameters is not yet supported'}, status=400)

    return JsonResponse([best_fit_ephemeris, cov], safe=False, status=200)
