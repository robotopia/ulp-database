from django.contrib import admin
from .models import *
from common.admin import PermissionFieldset

# Register your models here.

@admin.register(TimeOfArrival)
class TimeOfArrivalAdmin(admin.ModelAdmin):
    list_display = ['pk', 'ulp', 'mjd', 'telescope_name', 'raw_mjd', 'freq', 'barycentred', 'dedispersed']
    list_filter = ['ulp']
    fieldsets = [
        PermissionFieldset,
        (
            None, {
                "fields": ["ulp", ("mjd", "mjd_err"), "raw_mjd", "telescope_name", "freq", "bw", "spectral_index", "rotation_measure", "peak_flux_Jy", "upper_limit", "pulse_width", "barycentred", "dedispersed", "notes", "plots", "lightcurve", "pulse"],
            }
        ),
    ]

@admin.register(LightcurvePoint)
class LightcurvePointAdmin(admin.ModelAdmin):
    list_display = ['pk', 'lightcurve', 'sample_number', 'pol', 'value']
    list_filter = ['lightcurve__ulp']

@admin.register(Lightcurve)
class LightcurveAdmin(admin.ModelAdmin):
    list_display = ['pk', 'ulp', 'telescope', 't0', 't0_gps', 'freq']
    list_filter = ['ulp', 'telescope']
    fieldsets = [
        PermissionFieldset,
        (
            None, {
                "fields": ["ulp", "telescope", ("freq", "bw"), ("t0", "dt"), ("dm", "dm_freq")],
            }
        ),
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

@admin.register(WorkingEphemeris)
class WorkingEphemerisAdmin(admin.ModelAdmin):
    list_display = ['pk', 'owner', 'ulp', 'pepoch', 'p0', 'p1', 'pb', 'dm']
    list_filter = ['owner', 'ulp']
    fieldsets = [
        PermissionFieldset,
        (None, {"fields": ["ulp"]}),
        (
            "Ephemeris values", {
                "fields": ["pepoch", "p0", "p1", "pb", "dm"],
            }
        ),
        (
            "Spectrum", {
                "fields": ["spec_alpha", "spec_q"],
            }
        ),
    ]

@admin.register(Pulse)
class PulseAdmin(admin.ModelAdmin):
    list_display = ['pk', 'lightcurve', 'mjd_start', 'mjd_end', 'tags']
    list_filter = ['lightcurve__owner', 'lightcurve__ulp', 'lightcurve__telescope']

@admin.register(Template)
class TemplateAdmin(admin.ModelAdmin):
    list_display = ['pk', 'working_ephemeris', 'owner', 'updated']
    list_filter = ['working_ephemeris__ulp']
    fieldsets = [
        PermissionFieldset,
        (
            None, {
                "fields": ["working_ephemeris"],
            }
        ),
    ]

@admin.register(TemplateComponent)
class TemplateComponentAdmin(admin.ModelAdmin):
    list_display = ['pk', 'template', 'weight', 'mu', 'sigma']


@admin.register(Toa)
class ToaAdmin(admin.ModelAdmin):
    list_display = ['pk', 'toa_mjd', 'toa_err_s', 'template', 'pulse_number', 'ampl']
    list_filter = ['template__working_ephemeris__ulp']
