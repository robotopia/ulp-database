from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.contrib.auth import authenticate, login, logout
from django.db.models import Q, Max
from . import models
from decimal import Decimal

import pygedm
from astropy.coordinates import SkyCoord

# Create your views here.

def index(request):
    return render(request, 'published/index.html')

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

    parameter_set = get_object_or_404(models.ParameterSet, name='main_table')
    measurements = get_accessible_measurements(request, parameter_set=parameter_set)

    ulps = list({measurement.ulp for measurement in measurements})
    parameters = [parameter for parameter in parameter_set.parameters.all()]

    rows = {}
    for ulp in ulps:
        rows[ulp] = {}
        for parameter in parameters:
            rows[ulp][parameter.name] = measurements.filter(ulp=ulp, parameter=parameter).order_by('updated').last()

    context = {
        'rows': rows,
    }

    return render(request, 'published/ppdot.html', context)


def galactic_view(request):

    method = request.GET.get('method') or 'DM'
    model = request.GET.get('model') or 'ne2001'
    dm_dist_frac_err = float(request.GET.get('dm_dist_frac_err') or 0)

    parameter_set = get_object_or_404(models.ParameterSet, name="position")
    measurements = get_accessible_measurements(request, parameter_set=parameter_set)

    ulps = list({measurement.ulp for measurement in measurements})
    parameters = [parameter for parameter in parameter_set.parameters.all()]

    distance_parameter = models.Parameter.objects.get(name='Distance')

    # Organise measurements by ULP and parameter
    values = {}
    for ulp in ulps:
        values[ulp] = {}
        for parameter in parameters:
            values[ulp][parameter.name] = measurements.filter(ulp=ulp, parameter=parameter).order_by('updated').last()

        # Pull out some quantities for future convenience
        ra = values[ulp]['Right ascension'].astropy_quantity
        dec = values[ulp]['Declination'].astropy_quantity

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
        except:
            pass

        # Update the distance according to the method and the model
        if method == 'DM' and dm is not None:
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

    context = {
        'values': values,
        'method': method,
        'model': model,
        'dm_dist_frac_err': dm_dist_frac_err,
    }

    return render(request, 'published/galactic_view.html', context)

