from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404
from django.db.models import Q
from . import models

from published import models as published_models
from published.views import get_accessible_measurements

import numpy as np
import astropy.units as u

# Create your views here.

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


def timing_residual_view(request, pk):

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

    ephemeris_measurements = models.EphemerisMeasurement.objects.filter(
        measurement__owner=request.user,
        measurement__ulp=ulp,
    )

    if not ephemeris_measurements.exists():
        return HttpResponse(status=404)

    # Construct a dictionary out of the ephemeris
    ephemeris = {e.ephemeris_parameter.tempo_name: e.value for e in ephemeris_measurements}

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

    context = {
        'ulp': ulp,
        'ephemeris': ephemeris,
        'periods': periods,
    }

    mjds = [float(toa.mjd) for toa in toas]
    xdata_min = np.min(mjds)
    xdata_max = np.max(mjds)
    xdata_range = xdata_max - xdata_min
    plot_specs = {
        'xmin': xdata_min - 0.05*xdata_range,  # With an extra margin buffer
        'xmax': xdata_max + 0.05*xdata_range,
        'xrange': xdata_range,
    }
    context['plot_specs'] = plot_specs

    return render(request, 'data/timing_residuals.html', context)
