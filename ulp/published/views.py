from django.shortcuts import render
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from django.db.models import Q, Max
from . import models

# Create your views here.

def index(request):
    return HttpResponse("Hello, world. You're at the 'published' index.")

def main_table(request):

    parameter_set = get_object_or_404(models.ParameterSet, name='main_table')

    measurements = models.Measurement.objects.filter(parameter__in=parameter_set.parameters.all()).filter(
        Q(article__isnull=False) |  # It's published, and therefore automatically accessible by everyone
        Q(owner=request.user) |  # The owner can always see their own measurements
        Q(access=models.Measurement.ACCESS_PUBLIC) |  # Include measurements explicitly marked as public
        (Q(access=models.Measurement.ACCESS_GROUP) &  # But if it's marked as group-accessible...
         Q(access_groups__in=request.user.groups.all()))  # ...then the user must be in of the allowed groups.
    )

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

    print(context)

    return render(request, 'published/main_table.html', context)
