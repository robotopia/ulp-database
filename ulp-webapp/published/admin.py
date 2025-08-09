from django.contrib import admin
from .models import *

# Register your models here.

@admin.register(Article)
class ArticleAdmin(admin.ModelAdmin):
    list_display = ['citet_text', 'title', 'doi']

@admin.register(Covariance)
class CovarianceAdmin(admin.ModelAdmin):
    list_display = ['measurement1', 'measurement2', 'covariance']

@admin.register(FrequencyBand)
class FrequencyBandAdmin(admin.ModelAdmin):
    list_display = ['pk', 'name', 'symbol', 'range_display']

@admin.register(Measurement)
class MeasurementAdmin(admin.ModelAdmin):
    list_display = ['ulp', 'parameter', 'stokes', 'formatted_quantity', 'article', 'owner', 'accessible_by']
    list_filter = ['ulp', 'parameter', 'access', 'article']
    fieldsets = (
        ('Main', {'fields': ('ulp', 'parameter', ('quantity', 'power_of_10'), 'article', 'date',)}),
        ('Uncertainty', {'fields': ('err', ('err_lo', 'err_hi'), 'precision', ('error_sigma', 'error_sigma_type'), ('chisq', 'reduced_chisq')), 'classes': ('collapse',)}),
        ('Frequency', {'fields': ('freq_band', ('freq_lo', 'freq_ctr', 'freq_hi'), 'freq_astropy_units'), 'classes': ('collapse',)}),
        ('Polarisation', {'fields': ('stokes',), 'classes': ('collapse',)}),
        ('Publish', {'fields': ('owner', 'access', 'access_groups'), 'classes': ('collapse',)}),
        ('Display options', {'fields': (('lower_limit', 'upper_limit'), 'error_is_range', 'approximation', 'angle_display',), 'classes': ('collapse',)}),
        ('Other', {'fields': ('notes',)}),
    )

@admin.register(Parameter)
class ParameterAdmin(admin.ModelAdmin):
    list_display = ['name', 'ascii_symbol', 'astropy_unit']

@admin.register(ParameterSet)
class ParameterSetAdmin(admin.ModelAdmin):
    list_display = ['name']

@admin.register(Progenitor)
class ProgenitorAdmin(admin.ModelAdmin):
    list_display = ['pk', 'name', 'abbr']

@admin.register(ProgenitorClaim)
class ProgenitorClaimAdmin(admin.ModelAdmin):
    list_display = ['pk', 'ulp', '__str__', 'article']
    list_filter = ['ulp']

@admin.register(Ulp)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['pk', 'name', 'abbr', 'best_progenitor_claims']

@admin.register(UserSetting)
class UserSettingAdmin(admin.ModelAdmin):
    list_display = ['pk', 'user', 'site_theme']

