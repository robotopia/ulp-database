from django.shortcuts import render, redirect
from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404
from django.contrib.auth import authenticate, login, logout
from django.db.models import Q, Max
import requests
from . import models
from decimal import Decimal
import csv
import io
import json

import pygedm
from astropy.coordinates import SkyCoord

from psrcat import models as psrcat_models

# Create your views here.

def swap_index(request):
    print("Here")
    res = render(request, 'published/index.html')
    print("here2")
    return res

def index(request):
    return redirect('galactic_view')

def get_accessible_measurements(request, parameter_set=None, ulp=None):

    queryset = models.Measurement.objects.all()

    if parameter_set is not None:
        queryset = queryset.filter(parameter__in=parameter_set.parameters.all())

    if ulp is not None:
        queryset = queryset.filter(ulp=ulp)

    if request.user.is_authenticated:
        queryset = queryset.filter(
            Q(article__isnull=False) |  # It's published, and therefore automatically accessible by everyone
            Q(owner=request.user) |  # The owner can always see their own measurements
            Q(access=models.Measurement.ACCESS_PUBLIC) |  # Include measurements explicitly marked as public
            (Q(access=models.Measurement.ACCESS_GROUP) &  # But if it's marked as group-accessible...
             Q(access_groups__in=request.user.groups.all()))  # ...then the user must be in of the allowed groups.
        )
    else:
        queryset = queryset.filter(
            Q(article__isnull=False) |  # It's published, and therefore automatically accessible by everyone
            Q(access=models.Measurement.ACCESS_PUBLIC)  # Include measurements explicitly marked as public
        )

    return queryset


def parameter_set_table_view(request, pk):

    parameter_set = get_object_or_404(models.ParameterSet, pk=pk)
    measurements = get_accessible_measurements(request, parameter_set=parameter_set)

    ulps = list({measurement.ulp for measurement in measurements})
    parameters = [parameter for parameter in parameter_set.parameters.all()]

    rows = {}
    for ulp in ulps:
        rows[ulp] = {}
        for parameter in parameters:
            rows[ulp][parameter.name] = measurements.filter(ulp=ulp, parameter=parameter).order_by('updated').last()

    # Also get the list of parameter sets to populate the dropdown
    parameter_sets = models.ParameterSet.objects.all()

    context = {
        'parameters': parameters,
        'parameter_sets': parameter_sets,
        'rows': rows,
    }

    return render(request, 'published/main_table.html', context)


def ulp_view(request, pk):

    ulp = get_object_or_404(models.Ulp, pk=pk)

    measurements = get_accessible_measurements(request, ulp=ulp)

    if not measurements.exists():
        return HttpResponse(status=404)

    context = {
        'ulp': ulp,
        'measurements': measurements,
    }

    return render(request, 'published/ulp.html', context)


def ppdot_view(request):

    import_mcgill = request.GET.get("mcgill") == "on"

    context = {
        'import_mcgill': import_mcgill,
    }

    return render(request, 'published/ppdot.html', context)


def psrcat_data(request):

    pulsars = psrcat_models.Pulsar.objects.using('psrcat')
    validated_psrcat_json_data = []
    for pulsar in pulsars:
        F0s = pulsar.parameter_set.filter(parametertype__label='F', timederivative=0)
        F1s = pulsar.parameter_set.filter(parametertype__label='F', timederivative=1)
        if F0s.exists() and F1s.exists():
            P0 = 1/float(F0s.first().value) # TODO: Get the currently accepted value instead of the 'first' value
            P1 = -float(F1s.first().value)*P0**2 # TODO: Get the currently accepted value instead of the 'first' value
            validated_psrcat_json_data.append(
                {
                    'name': pulsar.jname,
                    'P': P0,
                    'Pdot': P1,
                    'Pdot_err': None,
                    'Pdot__upper_limit': False,
                }
            )

    return JsonResponse(validated_psrcat_json_data, safe=False)


def mcgill_data(request):

    url = "https://www.physics.mcgill.ca/~pulsar/magnetar/TabO1.csv"
    mcgill = requests.get(url, stream=True)
    reader = csv.DictReader(io.StringIO(mcgill.text))
    validated_mcgill_json_data = [
        {
            'name': magnetar['Name'],
            'P': float(magnetar['Period']) if magnetar['Period'] != "" else None,
            'Pdot': float(magnetar['Pdot']) if magnetar['Pdot'] != "" else None,
            'Pdot_err': float(magnetar['Pdot_Err']) if magnetar['Pdot_Err'] != "" else None,
            'Pdot__upper_limit': magnetar['Pdot_lim'] == '<',
            'radio': 'R' in magnetar['Bands'],
        }
        for magnetar in list(reader)
    ]

    return JsonResponse(validated_mcgill_json_data, safe=False)


def table_data(request, pk):

    parameter_set = get_object_or_404(models.ParameterSet, pk=pk)

    measurements = get_accessible_measurements(request, parameter_set=parameter_set)

    ulps = list({measurement.ulp for measurement in measurements})
    parameters = [parameter for parameter in parameter_set.parameters.all()]

    plot_data = []
    for ulp in ulps:
        ulp_dict = {'id': ulp.id, 'name': ulp.name}
        for parameter in parameters:
            latest_measurement = measurements.filter(ulp=ulp, parameter=parameter).order_by('updated').last()
            if latest_measurement is None:
                continue
            ulp_dict[parameter.ascii_symbol] = latest_measurement.astropy_quantity.value
            if latest_measurement.astropy_err is not None:
                ulp_dict[parameter.ascii_symbol + '_err'] = latest_measurement.astropy_err.value
            ulp_dict[parameter.ascii_symbol + '_unit'] = latest_measurement.parameter.astropy_unit
            ulp_dict[parameter.ascii_symbol + '__upper_limit'] = latest_measurement.upper_limit == True
            ulp_dict[parameter.ascii_symbol + '__lower_limit'] = latest_measurement.lower_limit == True
        ulp_dict['radio'] = True
        plot_data.append(ulp_dict)

    return JsonResponse(plot_data, safe=False)


def galactic_view(request):

    model = request.GET.get('model') or 'ne2001'
    dm_dist_frac_err = float(request.GET.get('dm_dist_frac_err') or 0.3)

    parameter_set = get_object_or_404(models.ParameterSet, name="position")
    measurements = get_accessible_measurements(request, parameter_set=parameter_set)

    ulps = list({measurement.ulp for measurement in measurements})
    parameters = [parameter for parameter in parameter_set.parameters.all()]

    distance_parameter = models.Parameter.objects.get(name='Distance')

    colours = {
        'published': 'white',
        'dm': 'yellow',
        'unknown': 'gray',
    }

    # Organise measurements by ULP and parameter
    values = {}
    for ulp in ulps:
        values[ulp] = {}
        for parameter in parameters:
            values[ulp][parameter.name] = measurements.filter(ulp=ulp, parameter=parameter).order_by('updated').last()

        # Pull out some quantities for future convenience
        try:
            ra = values[ulp]['Right ascension'].astropy_quantity
            dec = values[ulp]['Declination'].astropy_quantity
        except:
            # If they don't have a position, it cannot be plotted, so ignore this ULP
            values.pop(ulp)
            continue

        coord = SkyCoord(ra=ra, dec=dec, frame='icrs')
        values[ulp]['gal_long'] = coord.galactic.l.value
        values[ulp]['gal_lat'] = coord.galactic.b.value

        try:
            dm = values[ulp]['Dispersion measure'].astropy_quantity
            dm_err = values[ulp]['Dispersion measure'].astropy_err
        except:
            dm = None
            dm_err = None

        try:
            dist = values[ulp]['Distance'].astropy_quantity
            dist_err = values[ulp]['Distance'].astropy_err

            # Convert to galactic coordinates
            values[ulp]['dist'] = dist.to('kpc').value
            values[ulp]['near_dist'] = (dist - dist_err).to('kpc').value
            values[ulp]['far_dist'] = (dist + dist_err).to('kpc').value
            values[ulp]['colour'] = colours['published']
        except:
            # Update the distance according to the model
            if dm is not None:
                pygedm_dist, _ = pygedm.dm_to_dist(coord.galactic.l, coord.galactic.b, dm, method=model)
                if dm_dist_frac_err <= 0:
                    pygedm_dist_lo, _ = pygedm.dm_to_dist(coord.galactic.l, coord.galactic.b, dm - dm_err, method=model)
                    pygedm_dist_hi, _ = pygedm.dm_to_dist(coord.galactic.l, coord.galactic.b, dm + dm_err, method=model)
                else:
                    pygedm_dist_lo = pygedm_dist * (1 - dm_dist_frac_err/100) # dm_dist_frac_err are %'s
                    pygedm_dist_hi = pygedm_dist * (1 + dm_dist_frac_err/100) # dm_dist_frac_err are %'s
                values[ulp]['dist'] = pygedm_dist.to('kpc').value
                values[ulp]['near_dist'] = pygedm_dist_lo.to('kpc').value
                values[ulp]['far_dist'] = pygedm_dist_hi.to('kpc').value
                values[ulp]['colour'] = colours['dm']

        if values[ulp]['Distance'] is None and dm is None:
            values[ulp]['dist'] = 15
            values[ulp]['near_dist'] = 0
            values[ulp]['far_dist'] = 30
            values[ulp]['colour'] = colours['unknown']

    context = {
        'values': values,
        'model': model,
        'colours': colours,
        'dm_dist_frac_err': int(dm_dist_frac_err),
    }

    return render(request, 'published/galactic_view.html', context)

