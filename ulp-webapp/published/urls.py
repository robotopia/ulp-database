from django.urls import re_path, path, include

from . import views

urlpatterns = [
    re_path('^$', views.index, name="index"),
    re_path('^swap_index/$', views.swap_index, name="swap_index"),
    re_path(r'^tables/(?P<pk>[0-9]+)$', views.parameter_set_table_view, name="parameter_set_table_view"),
    re_path(r'^ulp/(?P<pk>[0-9]+)$', views.ulp_view, name="ulp"),
    re_path(r'^galactic_view$', views.galactic_view, name="galactic_view"),
    re_path(r'^ppdot$', views.ppdot_view, name="ppdot_view"),
    re_path(r'^user_settings$', views.user_settings, name="user_settings"),
    re_path(r'^api/table/(?P<pk>[0-9]+)', views.table_data, name="table_data"),
    re_path(r'^api/mcgill_table', views.mcgill_data, name="mcgill_data"),
    re_path(r'^api/psrcat_table', views.psrcat_data, name="psrcat_data"),
    path('', include('django.contrib.auth.urls')),
]
