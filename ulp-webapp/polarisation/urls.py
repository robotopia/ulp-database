from django.urls import re_path, path, include

from . import views
from common.utils import update_permissions

urlpatterns = [
    re_path(r'^$', views.index, name="index"),
]
