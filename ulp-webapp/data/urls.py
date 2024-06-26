from django.urls import re_path, path, include

from . import views

urlpatterns = [
    re_path(r'^timing/$', views.timing_choose_ulp_view, name="timing_choose_ulp"),
    re_path(r'^timing/(?P<pk>[0-9]+)$', views.timing_residual_view, name="timing_residuals"),
    re_path(r'^api/toa_data/(?P<pk>[0-9]+)', views.toa_data, name="toa_data"),
]
