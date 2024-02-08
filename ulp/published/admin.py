from django.contrib import admin
from .models import *

# Register your models here.

@admin.register(Ulp)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['name', 'abbr']

@admin.register(Article)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['citet_text', 'doi']

@admin.register(Parameter)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['name', 'ascii_symbol', 'unit']

@admin.register(Measurement)
class UlpAdmin(admin.ModelAdmin):
    list_display = ['ulp', 'article', 'parameter', 'formatted_quantity']

