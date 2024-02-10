from django.shortcuts import render
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.db.models import Q, Max
from . import models
from decimal import Decimal

import pygedm
from astropy.coordinates import SkyCoord

# Create your views here.

def index(request):
    return HttpResponse("Hello, world. You're at the 'published' index.")

def get_accessible_measurements(request, parameter_set):

    return models.Measurement.objects.filter(parameter__in=parameter_set.parameters.all()).filter(
        Q(article__isnull=False) |  # It's published, and therefore automatically accessible by everyone
        Q(owner=request.user) |  # The owner can always see their own measurements
        Q(access=models.Measurement.ACCESS_PUBLIC) |  # Include measurements explicitly marked as public
        (Q(access=models.Measurement.ACCESS_GROUP) &  # But if it's marked as group-accessible...
         Q(access_groups__in=request.user.groups.all()))  # ...then the user must be in of the allowed groups.
    )

def parameter_set_table_view(request, pk):

    parameter_set = get_object_or_404(models.ParameterSet, pk=pk)
    measurements = get_accessible_measurements(request, parameter_set)

    ulps = list({measurement.ulp for measurement in measurements})
    parameters = [parameter for parameter in parameter_set.parameters.all()]

    rows = {}
    for ulp in ulps:
        rows[ulp.name] = {}
        for parameter in parameters:
            rows[ulp.name][parameter.name] = measurements.filter(ulp=ulp, parameter=parameter).order_by('updated').last()

    context = {
        'parameters': parameters,
        'rows': rows,
    }

    return render(request, 'published/main_table.html', context)

def galactic_view(request):

    method = request.GET.get('method') or 'DM'
    model = request.GET.get('model') or 'ne2001'
    dm_dist_frac_err = float(request.GET.get('dm_dist_frac_err') or 0)

    parameter_set = get_object_or_404(models.ParameterSet, name="position")
    measurements = get_accessible_measurements(request, parameter_set)

    ulps = list({measurement.ulp for measurement in measurements})
    parameters = [parameter for parameter in parameter_set.parameters.all()]

    distance_parameter = models.Parameter.objects.get(name='Distance')

    # Organise measurements by ULP and parameter
    values = {}
    for ulp in ulps:
        values[ulp.name] = {}
        for parameter in parameters:
            values[ulp.name][parameter.name] = measurements.filter(ulp=ulp, parameter=parameter).order_by('updated').last()

        # Pull out some quantities for future convenience
        ra = values[ulp.name]['Right ascension'].astropy_quantity
        dec = values[ulp.name]['Declination'].astropy_quantity
        try:
            dm = values[ulp.name]['Dispersion measure'].astropy_quantity
            dm_err = values[ulp.name]['Dispersion measure'].astropy_err
        except:
            dm = None
            dm_err = None
        dist = values[ulp.name]['Distance'].astropy_quantity
        dist_err = values[ulp.name]['Distance'].astropy_err

        # Convert to galactic coordinates
        coord = SkyCoord(ra=ra, dec=dec, frame='icrs')
        values[ulp.name]['gal_long'] = coord.galactic.l.value
        values[ulp.name]['gal_lat'] = coord.galactic.b.value
        values[ulp.name]['dist'] = dist.to('kpc').value
        values[ulp.name]['near_dist'] = (dist - dist_err).to('kpc').value
        values[ulp.name]['far_dist'] = (dist + dist_err).to('kpc').value

        # Update the distance according to the method and the model
        if method == 'DM' and dm is not None:
            pygedm_dist, _ = pygedm.dm_to_dist(coord.galactic.l, coord.galactic.b, dm, method=model)
            print(coord.galactic.l, coord.galactic.b, dm, " -> ", pygedm_dist)
            if dm_dist_frac_err <= 0:
                pygedm_dist_lo, _ = pygedm.dm_to_dist(coord.galactic.l, coord.galactic.b, dm - dm_err, method=model)
                pygedm_dist_hi, _ = pygedm.dm_to_dist(coord.galactic.l, coord.galactic.b, dm + dm_err, method=model)
            else:
                pygedm_dist_lo = pygedm_dist * (1 - dm_dist_frac_err/100) # dm_dist_frac_err are %'s
                pygedm_dist_hi = pygedm_dist * (1 + dm_dist_frac_err/100) # dm_dist_frac_err are %'s
            values[ulp.name]['dist'] = pygedm_dist.to('kpc').value
            values[ulp.name]['near_dist'] = pygedm_dist_lo.to('kpc').value
            values[ulp.name]['far_dist'] = pygedm_dist_hi.to('kpc').value

    context = {
        'values': values,
        'method': method,
        'model': model,
        'dm_dist_frac_err': dm_dist_frac_err,
    }

    print(context)

    return render(request, 'published/galactic_view.html', context)
