from django.contrib import admin
from .models import *
from common.admin import PermissionFieldset

# Register your models here.

@admin.register(TimeOfArrival)
class TimeOfArrivalAdmin(admin.ModelAdmin):
    list_display = ['pk', 'ulp', 'mjd', 'telescope_name', 'raw_mjd', 'freq', 'barycentred', 'dedispersed']
    list_filter = ['ulp']
    fieldsets = [
        (
            None, {
                "fields": ["ulp", ("mjd", "mjd_err"), "raw_mjd", "telescope_name", "freq", "bw", "spectral_index", "rotation_measure", "peak_flux_Jy", "upper_limit", "pulse_width", "barycentred", "dedispersed", "notes", "plots",],
            }
        ),
        PermissionFieldset,
    ]

@admin.register(EphemerisParameter)
class EphemerisParameterAdmin(admin.ModelAdmin):
    list_display = ['tempo_name', 'parameter', 'tempo_astropy_units']

@admin.register(EphemerisMeasurement)
class EphemerisMeasurementAdmin(admin.ModelAdmin):
    list_display = ['pk', 'ulp', 'ephemeris_parameter', 'owner', 'value', 'measurement',]
    list_filter = ['measurement__ulp', 'ephemeris_parameter', 'measurement__owner']

@admin.register(Plot)
class PlotAdmin(admin.ModelAdmin):
    list_display = ['pk', 'image', 'owner']
    list_filter = ['owner']

