from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="index"),
    path("main_table", views.main_table, name="main_table"),
]
