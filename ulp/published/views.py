from django.shortcuts import render
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from . import models

# Create your views here.

def index(request):
    return HttpResponse("Hello, world. You're at the 'published' index.")

def main_table(request):

    parameter_set = get_object_or_404(models.ParameterSet, name='main_table')
    context = {
        'parameter_set': parameter_set,
    }

    return render(request, 'published/main_table.html', context)
