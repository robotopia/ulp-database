from django.urls import re_path, path, include

from . import views

urlpatterns = [
    re_path(r'^timing/$', views.timing_choose_ulp_view, name="timing_choose_ulp"),
    re_path(r'^timing/(?P<pk>[0-9]+)$', views.timing_residual_view, name="timing_residuals"),
    re_path(r'^timing/toa/(?P<pk>[0-9]+)$', views.toa_detail_view, name="toa_detail_view"),
    re_path(r'^timing/toas/(?P<pk>[0-9]+)$', views.toas_view, name="toas_view"),
    re_path(r'^api/toa_data/(?P<pk>[0-9]+)', views.toa_data, name="toa_data"),
    re_path(r'^api/update_toa$', views.update_toa, name="update_toa"),
]
