from django.db.models import Q
from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.apps import apps
from django.contrib.auth.models import User, Group
import json
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import data.models as data_models

def permitted_to_view_filter(queryset, user):

    return queryset.filter(
        Q(owner=user) |
        Q(can_edit_groups__user=user) |
        Q(can_edit_users=user) |
        Q(can_view_groups__user=user) |
        Q(can_view_users=user)
    ).distinct()


def permitted_to_edit_filter(queryset, user):

    return queryset.filter(
        Q(owner=user) |
        Q(can_edit_groups__user=user) |
        Q(can_edit_users=user)
    ).distinct()


def permitted_to_delete_filter(queryset, user):

    return queryset.filter(Q(owner=user)).distinct()


def update_permissions(request):

    # First of all, they have to be logged in
    if not request.user.is_authenticated:
        return HttpResponse(status=404)

    # Turn the data into a dictionary
    try:
        data = json.loads(request.body.decode('utf-8'))
        app_name = data['app']
        model_name = data['model']
        pk = int(data['pk'])
        group_or_user = data['group_or_user']
        name = data['name']
        permission_type = data['permission_type']
        permit = bool(data['permit'])

        model = apps.get_model(app_name, model_name)
        obj = get_object_or_404(model, pk=pk)

    except:
        return HttpResponse(status=400)

    if not request.user == obj.owner:
        return HttpResponse(status=401)

    if group_or_user == 'group':
        group = get_object_or_404(Group, name=name)
        if permission_type == 'view':
            if permit:
                obj.can_view_groups.add(group)
            else:
                obj.can_view_groups.remove(group)
        elif permission_type == 'edit':
            if permit:
                obj.can_edit_groups.add(group)
            else:
                obj.can_edit_groups.remove(group)
    elif group_or_user == 'user':
        user = get_object_or_404(User, username=username)
        if permission_type == 'view':
            if permit:
                obj.can_view_users.add(user)
            else:
                obj.can_view_users.remove(user)
        elif permission_type == 'edit':
            if permit:
                obj.can_edit_users.add(user)
            else:
                obj.can_edit_users.remove(user)

    return HttpResponse(status=200)


def barycentre(ulp, times, location):
    '''
    ulp is a Ulp object
    times is an array-like object of MJDs (can be Time array)
    '''

    ra, dec = [
        data_models.EphemerisMeasurement.objects.filter(
            ephemeris_parameter__tempo_name=param,
            measurement__ulp=ulp
            ).first().measurement.quantity for param in ['RAJ', 'DECJ']
    ]

    direction = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    location_times = Time(times, format='mjd', scale='utc', location=location) # Tag times with location
    bc_correction = location_times.light_travel_time(direction, ephemeris='jpl') # Calculate correction
    bc_times = (location_times.tdb + bc_correction).value # Apply correction, and convert back to MJD values
    return bc_times
