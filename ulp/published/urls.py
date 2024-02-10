from django.urls import re_path

from . import views

urlpatterns = [
    re_path('^$', views.index, name="index"),
    re_path(r'^tables/(?P<pk>[0-9]+)$', views.parameter_set_table_view, name="parameter_set_table_view"),
    re_path(r'^galactic_view$', views.galactic_view, name="galactic_view"),
]
