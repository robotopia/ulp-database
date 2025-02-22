from django.urls import re_path, path, include

from . import views
from common.utils import update_permissions

urlpatterns = [
    re_path(r'^timing/$', views.timing_choose_ulp_view, name="timing_choose_ulp"),
    re_path(r'^timing/(?P<pk>[0-9]+)$', views.timing_residual_view, name="timing_residuals"),
    re_path(r'^timing/toa_detail/(?P<pk>[0-9]+)$', views.toa_detail_view, name="toa_detail_view"),
    re_path(r'^timing/toas/(?P<pk>[0-9]+)$', views.toas_view, name="toas_view"),
    re_path(r'^timing/toa/(?P<pk>[0-9]+)$', views.toa_view, name="toa_view"),  # Different sort of ToA to the above ones!!
    re_path(r'^timing/refit_toa/(?P<pk>[0-9]+)$', views.refit_toa, name="refit_toa"),
    re_path(r'^timing/lightcurve/(?P<pk>[0-9]+)$', views.lightcurve_view, name="lightcurve_view"),
    re_path(r'^timing/lightcurve_add/(?P<pk>[0-9]+)$', views.lightcurve_add, name="lightcurve_add"),
    re_path(r'^timing/folding/(?P<pk>[0-9]+)$', views.folding_view, name="folding_view"),
    re_path(r'^timing/folding_toa/(?P<pk>[0-9]+)$', views.folding_toa_view, name="folding_toa_view"),
    re_path(r'^timing/pulsestack/(?P<pk>[0-9]+)$', views.pulsestack_view, name="pulsestack_view"),
    re_path(r'^timing/update_working_ephemeris/(?P<pk>[0-9]+)$', views.update_working_ephemeris, name="update_working_ephemeris"),
    re_path(r'^timing/add_or_update_pulse/(?P<pk>[0-9]+)$', views.add_or_update_pulse, name="add_or_update_pulse"),
    re_path(r'^timing/toa_for_pulse/(?P<pk>[0-9]+)$', views.toa_for_pulse, name="toa_for_pulse"),
    re_path(r'^timing/download_toas/(?P<pk>[0-9]+)$', views.download_toas, name="download_toas"),
    re_path(r'^timing/download_working_ephemeris/(?P<pk>[0-9]+)$', views.download_working_ephemeris, name="download_working_ephemeris"),
    re_path(r'^api/toa_data/(?P<pk>[0-9]+)', views.toa_data, name="toa_data"),
    re_path(r'^api/update_toa$', views.update_toa, name="update_toa"),
    re_path(r'^api/update_permissions$', views.update_permissions, name="update_permissions"),
    re_path(r'^api/upload_lightcurve$', views.upload_lightcurve, name="upload_lightcurve"),
]
