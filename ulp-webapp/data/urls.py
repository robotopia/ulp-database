from django.urls import re_path, path, include

from . import views

urlpatterns = [
    #re_path('^$', views.index, name="index"),
    #re_path(r'^tables/(?P<pk>[0-9]+)$', views.parameter_set_table_view, name="parameter_set_table_view"),
    re_path(r'^timing/(?P<pk>[0-9]+)$', views.timing_residual_view, name="timing_residuals"),
    #re_path(r'^galactic_view$', views.galactic_view, name="galactic_view"),
    #re_path(r'^ppdot$', views.ppdot_view, name="ppdot_view"),
    #re_path(r'^api/table/(?P<pk>[0-9]+)', views.table_data, name="table_data"),
    #re_path(r'^api/mcgill_table', views.mcgill_data, name="mcgill_data"),
    #path('', include('django.contrib.auth.urls')),
]
