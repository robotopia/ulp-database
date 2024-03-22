from django.shortcuts import render
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.db.models import Q
from . import models

from published import models as published_models
from published.views import get_accessible_measurements

import numpy as np
import astropy.units as u

# Create your views here.

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
        'periods': periods,
    }

    # The folding period is obtained from the request body
    if request.method == 'POST':
        folding_period = request.POST.get('folding_period')
    else:
        folding_period = periods.first()

    # It's still possible that there is no selected folding period, in which
    # case no plot will be shown. But if there is a folding period, calculate
    # the residuals

    if folding_period is not None:

        # For the moment, choose an arbitrary PEPOCH
        pepoch = 60000.0*u.day

        toas_mjd = np.array([toa.mjd for toa in toas]) * u.day
        toas_mjd_err = np.array([toa.mjd_err for toa in toas]) * u.day

        phases = ((toas_mjd - pepoch) / folding_period.astropy_quantity).decompose()
        # Calculate the residuals
        print('here01')
        residuals_percent = (phases*100 + 50) % 100 - 50
        # Subtract the mean
        print('here02')
        residuals_percent -= np.mean(residuals_percent)
        print('here03')

        # Associate the residuals with their corresponding TOAs
        plot_data = [
            {
                'x': (toas_mjd[i] - pepoch).value,
                'y': residuals_percent[i].value,
                'yerr': (toas.mjd_err[i]*u.day / folding_period.astropy_quantity).decompose().value,
            }
            for i in range(len(toas_mjd))
        ]

        # Get some bounds for the plot. The SVG viewBox is defined in terms of
        #   MJD since PEPOCH for the x-axis
        #   Residual as a % for the y-axis
        xs = [plot_data[i]['x'] for i in range(len(plot_data))]
        ys = [plot_data[i]['y'] for i in range(len(plot_data))]
        xdata_min = np.min(xs)
        xdata_max = np.max(xs)
        ydata_min = np.min(ys)
        ydata_max = np.max(ys)
        xdata_range = xdata_max - xdata_min
        ydata_range = ydata_max - ydata_min
        plot_specs = {
            'xmin': xdata_min - 0.05*xdata_range,  # With an extra margin buffer
            'ymin': ydata_min - 0.05*ydata_range,
            'xmax': xdata_max + 0.05*xdata_range,
            'ymax': ydata_max + 0.05*ydata_range,
            'xrange': xdata_range,
            'yrange': ydata_range,
        }
        print(plot_specs)

        # Add plot info to the context
        context['plot_specs'] = plot_specs
        context['folding_period'] = folding_period
        context['plot_data'] = plot_data
        context['pepoch'] = pepoch.value

    return render(request, 'data/timing_residuals.html', context)
