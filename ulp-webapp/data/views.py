from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404
from django.db.models import Q
from . import models

from published import models as published_models
from published.views import get_accessible_measurements

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz, get_sun
from astropy.constants import c
from astropy.time import Time

from spiceypy.spiceypy import spkezr, furnsh, j2000, spd, unload

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


def toa_data(request, pk):
    # Retrieve the selected ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    # Make sure the user has the permissions to view this ULP

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=404)

    # Second, they have to belong to a group that has been granted access to
    # this ULP's data
    if not ulp.data_access_groups.filter(user=request.user).exists():
        return HttpResponse(status=404)

    # Otherwise, grant them access, and get the TOAs!
    toas = models.TimeOfArrival.objects.filter(ulp=ulp)
    if not toas.exists():
        return HttpResponse(status=404)

    toas_json = [
        {
            'mjd': float(toa.mjd),
            'mjd_err': float(toa.mjd_err),
        } for toa in toas
    ]

    return JsonResponse(toas_json, safe=False)


def timing_choose_ulp_view(request):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=404)

    # Get ULPs to which they have access
    ulps = published_models.Ulp.objects.filter(
        data_access_groups__user=request.user
    )

    context = {
        'ulps': ulps,
    }

    return render(request, 'data/timing_choose_ulp.html', context)


def timing_residual_view(request, pk):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=404)

    # Retrieve the selected ULP
    ulp = get_object_or_404(published_models.Ulp, pk=pk)

    context = {'ulp': ulp}

    # Make sure the user has the permissions to view this ULP

    # Second, they have to belong to a group that has been granted access to
    # this ULP's data
    if not ulp.data_access_groups.filter(user=request.user).exists():
        return HttpResponse(status=404)

    # Otherwise, grant them access, and get the TOAs!
    toas = models.TimeOfArrival.objects.filter(ulp=ulp)
    if not toas.exists():
        return HttpResponse(status=404)

    ephemeris_measurements = models.EphemerisMeasurement.objects.filter(
        measurement__access_groups__user=request.user,
        measurement__ulp=ulp,
    )

    if not ephemeris_measurements.exists():
        return HttpResponse(status=404)

    # Construct a dictionary out of the ephemeris
    ephemeris = {e.ephemeris_parameter.tempo_name: e.value for e in ephemeris_measurements}

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
        Q(access=published_models.Measurement.ACCESS_PUBLIC) |  # Include measurements explicitly marked as public
        (Q(access=published_models.Measurement.ACCESS_GROUP) &  # But if it's marked as group-accessible...
         Q(access_groups__in=request.user.groups.all()))  # ...then the user must be in of the allowed groups.
    )
    context['periods'] = periods

    # Calculate some sensible initial plot dimensions
    mjds = Time([float(toa.mjd) for toa in toas], format='mjd')
    mjd_errs = [float(toa.mjd_err) for toa in toas] * u.d

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
    context['telescopes'] = EarthLocation.get_site_names()
    context['min_el'] = min_el
    context['max_sun_el'] = max_sun_el

    # For display purposes, convert RAJ and DECJ to hexagesimal format
    ephemeris['RAJ'] = Angle(ephemeris['RAJ'], unit=u.deg).to_string(unit=u.hour, sep=':')
    ephemeris['DECJ'] = Angle(ephemeris['DECJ'], unit=u.deg).to_string(sep=':')

    return render(request, 'data/timing_residuals.html', context)
