from django.shortcuts import render
from django.http import HttpResponse
from django.shortcuts import get_object_or_404

from published import models as published_models
from published.views import get_accessible_measurements

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

    context = {
        'ulp': ulp,
        'toas': toas,
    }

    return render(request, 'data/timing_residuals.html', context)
