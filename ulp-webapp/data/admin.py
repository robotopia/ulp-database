from django.contrib import admin
from .models import *

# Register your models here.

@admin.register(Telescope)
class TelescopeAdmin(admin.ModelAdmin):
    list_display = ['name', 'abbr']

@admin.register(Backend)
class BackendAdmin(admin.ModelAdmin):
    list_display = ['name', 'telescope', 'freq_ctr', 'bw', 'freq_units']

@admin.register(TimeOfArrival)
class TimeOfArrivalAdmin(admin.ModelAdmin):
    list_display = ['ulp', 'mjd', 'mjd_err', 'barycentred', 'dedispersed']
    list_filter = ['ulp']

@admin.register(EphemerisParameter)
class EphemerisParameterAdmin(admin.ModelAdmin):
    list_display = ['tempo_name', 'parameter', 'tempo_astropy_units']

@admin.register(EphemerisMeasurement)
class EphemerisMeasurementAdmin(admin.ModelAdmin):
    list_display = ['pk', 'ulp', 'ephemeris_parameter', 'owner', 'value', 'measurement',]
    list_filter = ['measurement__ulp', 'ephemeris_parameter', 'measurement__owner']

