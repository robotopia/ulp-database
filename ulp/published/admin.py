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
    list_display = ['ulp', 'parameter', 'formatted_quantity', 'access', 'article', 'updated']
    list_filter = ['ulp', 'parameter', 'access', 'article']

@admin.register(Parameter)
class ParameterAdmin(admin.ModelAdmin):
    list_display = ['name', 'ascii_symbol', 'astropy_unit']

@admin.register(ParameterSet)
class ParameterSetAdmin(admin.ModelAdmin):
    list_display = ['name']

@admin.register(Ulp)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['name', 'abbr']

