from django.shortcuts import render, get_object_or_404, redirect
from django.http import HttpResponse, JsonResponse
from django.db.models import Q, Value, BooleanField
from django.urls import reverse
from django.core.exceptions import ValidationError
from . import models
from common.utils import *

from published import models as published_models
from published.views import get_accessible_measurements

import numpy as np
import json
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz, get_sun
from astropy.constants import c
from astropy.time import Time
from decimal import Decimal

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

def calc_dmdelay(dm, flo, fhi):
    return 4.148808e3 * u.s * (dm.to('pc/cm3').value) * (1/flo.to('MHz').value**2 - 1/fhi.to('MHz').value**2)

def ephemeris_to_skycoord(ephemeris):
    '''
    ephemeris_measurements_qs is a QuerySet. If it contains both RAJ and DECJ,
    a SkyCoord object will be constructed therefrom.
    '''
    try:
        coord = SkyCoord(ra=ephemeris['RAJ'], dec=ephemeris['DECJ'], unit=(u.deg, u.deg), frame='icrs')
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

    return delays.to('s')


def calc_pulse_phase(time, ephemeris):
    '''
    time is an astropy Time object
    ephemeris['PEPOCH'] should be in days
    ephemeris['P0'] should be in seconds
    '''
    pepoch = Time(ephemeris['PEPOCH'], format='mjd')
    P0 = ephemeris['P0']*u.s
    return ((time - pepoch)/P0).decompose()


def calc_mjd(pulse_phase, ephemeris):
    '''
    This is the inverse of calc_pulse_phase()
    '''
    pepoch = Time(ephemeris['PEPOCH'], format='mjd')
    P0 = ephemeris['P0']*u.s
    return pepoch + pulse_phase*P0


def generate_toas(time_start, time_end, ephemeris):

    pulse_phase_start = calc_pulse_phase(time_start, ephemeris)
    pulse_phase_end = calc_pulse_phase(time_end, ephemeris)

    pulse_phases = np.arange(np.ceil(pulse_phase_start), pulse_phase_end)
    mjds = calc_mjd(pulse_phases, ephemeris)

    return mjds


def barycentre_toas(toas, coord):
    '''
    toa is a queryset containing toa objects

    This will check if the toa has been barycentred, and, if not,
    it will barycentre it.
    '''

    mjds = Time([float(toa.mjd) for toa in toas], format='mjd')
    corrections = bc_corr(coord, mjds)
    for i in range(len(toas)):
        toa = toas[i]
        if toa.raw_mjd is not None and not toa.barycentred:
            toa.mjd = toa.raw_mjd + Decimal(corrections[i].to('day').value)
            toa.barycentred = True
            toa.save()


def toa_data(request, pk):
    # Retrieve the selected ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Make sure the user has the permissions to view this ULP

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

    # Second, they have to belong to a group that has been granted access to
    # this ULP's data
    if not ulp.data_access_groups.filter(user=request.user).exists() and not request.user in ulp.whitelist_users.all():
        return HttpResponse(status=404)

    # Otherwise, grant them access, and get the TOAs!
    toas = models.TimeOfArrival.objects.filter(ulp=ulp)
    if not toas.exists():
        return HttpResponse(status=404)

    toas_json = [
        {
            'mjd': float(toa.mjd),
            'mjd_err': float(toa.mjd_err),
            'detail_link': reverse('toa_detail_view', args=[toa.pk]),
        } for toa in toas
    ]

    return JsonResponse(toas_json, safe=False)


def timing_choose_ulp_view(request):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

    # Get ULPs to which they have access
    ulps = published_models.Ulp.objects.filter(
        Q(data_access_groups__user=request.user) |
        Q(whitelist_users=request.user)
    ).distinct()

    context = {
        'ulps': ulps,
    }

    return render(request, 'data/timing_choose_ulp.html', context)


def timing_residual_view(request, pk):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

    # Retrieve the selected ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    context = {'ulp': ulp}

    # Make sure the user has the permissions to view this ULP

    # Second, they have to belong to a group that has been granted access to
    # this ULP's data
    if not ulp.data_access_groups.filter(user=request.user).exists() and not request.user in ulp.whitelist_users.all():
        return HttpResponse(status=404)

    # Otherwise, grant them access, and get the TOAs!
    toas = models.TimeOfArrival.objects.filter(ulp=ulp)
    if not toas.exists():
        return HttpResponse(status=404)

    ephemeris_measurements = models.EphemerisMeasurement.objects.filter(
        Q(measurement__ulp=ulp) &
        (Q(measurement__access_groups__user=request.user) |
         Q(measurement__ulp__whitelist_users=request.user)),
    )

    if not ephemeris_measurements.exists():
        return HttpResponse(status=404)

    # Construct a dictionary out of the ephemeris
    ephemeris = {e.ephemeris_parameter.tempo_name: e.value for e in ephemeris_measurements}

    ### WARNING: This is commented out deliberately. Only uncomment if you need to "manually"
    ### force a bunch of topocentric TOAs to be made barycentric. This makes changes
    ### to the database itself.
    coord = ephemeris_to_skycoord(ephemeris)
    barycentre_toas(toas, coord)

    # Do something similar for toas that are not yet dediserpsed, but for which a frequency is given
    if 'DM' in ephemeris.keys():
        dm = ephemeris['DM'] * u.pc / u.cm**3
        for toa in toas:
            if toa.freq is not None and toa.dedispersed == False:
                delay = calc_dmdelay(dm, toa.freq*u.MHz, np.inf*u.MHz)
                toa.mjd -= Decimal(delay.to('d').value)
                toa.dedispersed = True
                toa.save()

    output_toa_format = 'mjd'
    min_el = 0.0
    max_sun_el = 90.0

    if request.method == "POST":

        # Get form values
        PEPOCH = float(request.POST.get('pepoch', '0.0'))
        P0 = float(request.POST.get('folding-period', '0.0'))
        DM = float(request.POST.get('dm', '0.0'))
        RAJ = Angle(request.POST.get('raj', '00:00:00') + ' h').deg
        DECJ = Angle(request.POST.get('decj', '00:00:00') + ' d').deg

        mjd_start = Time(request.POST.get('mjd-start'), format='isot')
        mjd_end = Time(request.POST.get('mjd-end'), format='isot')
        mjd_range = Time([request.POST.get('mjd-start'), request.POST.get('mjd-end')], format='isot')

        context['mjd_start'] = mjd_start.isot
        context['mjd_end'] = mjd_end.isot

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
        ephemeris['PEPOCH'] = PEPOCH
        ephemeris['P0'] = P0
        ephemeris['DM'] = DM
        ephemeris['RAJ'] = RAJ
        ephemeris['DECJ'] = DECJ

        # If they've also provided other form values, make a table of predicted values
        if mjd_start is not None and mjd_end is not None and mjd_dispersion_frequency is not None and PEPOCH is not None and P0 is not None and telescope is not None:
            # First, assume the given mjd_start and mjd_end are in fact topocentric dispersed MJDs,
            # so to get the right range, convert them to dedispersed, barycentric
            coord = ephemeris_to_skycoord(ephemeris)
            dmdelay = calc_dmdelay(ephemeris['DM']*u.pc/u.cm**3, mjd_dispersion_frequency*u.MHz, np.inf*u.MHz)
            mjd_range -= dmdelay
            mjd_range += bc_corr(coord, mjd_range)

            predicted_barycentric_toas = generate_toas(mjd_range[0], mjd_range[1], ephemeris)
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

            context['predicted_toas'] = predicted_toas

        context['mjd_dispersion_frequency'] = mjd_dispersion_frequency

    context['ephemeris'] = ephemeris

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

    mjds = Time([float(toa.mjd) for toa in toas], format='mjd')
    mjd_errs = [float(toa.mjd_err) for toa in toas] * u.d

    # Calculate some sensible initial plot dimensions
    xdata_min = np.min(mjds).mjd
    xdata_max = np.max(mjds).mjd
    xdata_range = xdata_max - xdata_min

    min_phases = (calc_pulse_phase(mjds - mjd_errs, ephemeris) + 0.5) % 1 - 0.5
    max_phases = (calc_pulse_phase(mjds + mjd_errs, ephemeris) + 0.5) % 1 - 0.5
    ydata_min = np.min(min_phases)
    ydata_max = np.max(max_phases)
    ydata_range = ydata_max - ydata_min

    plot_specs = {
        'xmin': xdata_min - 0.05*xdata_range,  # With an extra margin buffer
        'xmax': xdata_max + 0.05*xdata_range,
        'xrange': xdata_range,
        'ymin': ydata_min - 0.05*ydata_range,
        'ymax': ydata_max + 0.05*ydata_range,
        'yrange': ydata_range,
    }
    context['plot_specs'] = plot_specs

    # Generate list of output TOA formats:
    context['output_toa_formats'] = list(Time.FORMATS)
    context['selected_output_toa_format'] = output_toa_format
    context['telescopes'] = site_names
    context['min_el'] = min_el
    context['max_sun_el'] = max_sun_el

    # For display purposes, convert RAJ and DECJ to hexagesimal format
    ephemeris['RAJ'] = Angle(ephemeris['RAJ'], unit=u.deg).to_string(unit=u.hour, sep=':')
    ephemeris['DECJ'] = Angle(ephemeris['DECJ'], unit=u.deg).to_string(sep=':')

    return render(request, 'data/timing_residuals.html', context)


def toa_detail_view(request, pk):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

    # Retrieve the selected ULP
    toa = get_object_or_404(models.TimeOfArrival, pk=pk)

    context = {
        'toa': toa,
        'freq_units': toa_freq_units,
    }

    return render(request, 'data/toa_detail.html', context)


def toas_view(request, pk):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

    # Retrieve the selected ToAs from the specified ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)
    toas = ulp.times_of_arrival.all()

    # Now limit only to those that this user has view access to
    toas = permitted_to_view_filter(toas, request.user)

    # Make sure the requested display units are dimensionally correct; if not, use default
    try:
        freq_units = request.GET.get("freq_units")
        foo = (1*u.Unit(freq_units)).to(toa_freq_units)
    except:
        freq_units = "MHz"

    try:
        mjd_err_units = request.GET.get("mjd_err_units") or "min"
        foo = (1*u.Unit(mjd_err_units)).to(toa_mjd_err_units)
    except:
        mjd_err_units = "min"

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


def add_or_update_pulse(request, pk):

    lightcurve = get_object_or_404(models.Lightcurve, pk=pk)

    if request.method == "POST":

        # I should be validating data here...

        if 'pulse_id' in request.POST.keys():
            # We're updating an existing pulse
            pulse = get_object_or_404(models.Pulse, pk=int(request.POST['pulse_id']))

            if request.POST['action'] == "Save":
                pulse.mjd_start = request.POST['mjd_start']
                pulse.mjd_end = request.POST['mjd_end']
                pulse.tags = request.POST['tags']
                pulse.save()
            elif request.POST['action'] == "Delete":
                print(f"about to delete {pulse}...")
                pulse.delete()

        else:
            # We're creating a new pulse
            pulse = models.Pulse(
                lightcurve=lightcurve,
                mjd_start=request.POST['mjd_start'],
                mjd_end=request.POST['mjd_end'],
                tags=request.POST['tags'],
            )

            pulse.save()

    return redirect('lightcurve_view', pk=pk)

def update_toa(request):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

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


def lightcurve_view(request, pk):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

    # Get the relevant Lightcurve object
    lightcurve = get_object_or_404(models.Lightcurve, pk=pk)

    # The user must be allowed to view this object
    if not lightcurve.can_view(request.user):
        return HttpResponse(status=403)

    context = {
        'lightcurve': lightcurve,
    }

    lightcurve_points = lightcurve.points.all()

    if lightcurve_points.exists():

        # An attempt to make the plot client-side...
        pols = list({p.pol for p in lightcurve.points.all()})
        data = [
            {
                "x": list(lightcurve.times(pol)),
                "y": list(lightcurve.values(pol)),
                "name": pol,
            } for pol in pols
        ]
        context['data'] = data

    return render(request, 'data/lightcurve.html', context)


def lightcurve_add(request, pk):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

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

        toas = [get_object_or_404(models.TimeOfArrival, pk=toa_pk) for toa_pk in request.POST.getlist('toas')]

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

        # Associate the TOAs
        for toa in toas:
            toa.lightcurve = lightcurve
            toa.save()

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


def folding_view(request, pk):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=401)

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
            'mjds': list(lc.bary_times()),
            'values': list(values + i), # The +i is for stacking
            'link': reverse('lightcurve_view', args=[lc.pk]),
            'pulses': [list(barycentre(ulp, [p.mjd_start, p.mjd_end], EarthLocation.of_site(lc.telescope))) for p in lc.pulses.all()],
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

    # If this is a POST request, save the provided values
    if request.method == "POST":
        for field in ['pepoch', 'p0', 'p1', 'pb', 'dm']:
            try:
                setattr(working_ephemeris, field, float(request.POST[field]))
            except:
                pass

        working_ephemeris.save()

    return render(request, 'data/folding.html', context)

