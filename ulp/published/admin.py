from django.contrib import admin
from .models import *

# Register your models here.

@admin.register(Article)
class ArticleAdmin(admin.ModelAdmin):
    list_display = ['citet_text', 'title', 'doi']

@admin.register(Covariance)
class CovarianceAdmin(admin.ModelAdmin):
    list_display = ['measurement1', 'measurement2', 'covariance']
    #list_filter = ['measurement

@admin.register(Measurement)
class MeasurementAdmin(admin.ModelAdmin):
    list_display = ['ulp', 'parameter', 'formatted_quantity', 'article', 'updated']
    list_filter = ['ulp', 'parameter', 'access', 'article']
    fieldsets = (
        ('Main', {'fields': ('ulp', 'parameter', ('quantity', 'power_of_10'), 'article', 'date',)}),
        ('Uncertainty', {'fields': ('err', ('err_lo', 'err_hi'), 'error_sigma', ('chisq', 'reduced_chisq')), 'classes': ('collapse',)}),
        ('Frequency', {'fields': (('freq_lo', 'freq_ctr', 'freq_hi'), 'freq_astropy_units'), 'classes': ('collapse',)}),
        ('Publish', {'fields': ('owner', 'access', 'access_groups'), 'classes': ('collapse',)}),
        ('Display options', {'fields': ('precision', ('lower_limit', 'upper_limit'), 'error_is_range', 'approximation', 'angle_display',), 'classes': ('collapse',)}),
        ('Other', {'fields': ('notes',)}),
    )

@admin.register(Parameter)
class ParameterAdmin(admin.ModelAdmin):
    list_display = ['name', 'ascii_symbol', 'astropy_unit']

@admin.register(ParameterSet)
class ParameterSetAdmin(admin.ModelAdmin):
    list_display = ['name']

@admin.register(Ulp)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['name', 'abbr']

