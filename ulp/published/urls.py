from django.urls import re_path, path, include

from . import views

urlpatterns = [
    re_path('^$', views.index, name="index"),
    re_path(r'^tables/(?P<pk>[0-9]+)$', views.parameter_set_table_view, name="parameter_set_table_view"),
    re_path(r'^ulp/(?P<pk>[0-9]+)$', views.ulp_view, name="ulp"),
    re_path(r'^galactic_view$', views.galactic_view, name="galactic_view"),
    path('', include('django.contrib.auth.urls')),
]
