from django.contrib import admin
from .models import *

# Register your models here.

@admin.register(Article)
class ArticleAdmin(admin.ModelAdmin):
    list_display = ['citet_text', 'doi']

@admin.register(Measurement)
class MeasurementAdmin(admin.ModelAdmin):
    list_display = ['ulp', 'article', 'parameter', 'formatted_quantity', 'access', 'updated']

@admin.register(Parameter)
class ParameterAdmin(admin.ModelAdmin):
    list_display = ['name', 'ascii_symbol', 'astropy_unit']

@admin.register(ParameterSet)
class ParameterSetAdmin(admin.ModelAdmin):
    list_display = ['name']

@admin.register(Ulp)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['name', 'abbr']

