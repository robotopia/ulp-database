from django.db.models import Q

def permitted_to_view(queryset, user):

    return queryset.filter(
        Q(owner=user) |
        Q(can_edit_groups__users__user=user) |
        Q(can_edit_users__user=user) |
        Q(can_view_groups__users__user=user) |
        Q(can_view_users__user=user)
    ).distinct()


def permitted_to_edit(queryset, user):

    return queryset.filter(
        Q(owner=user) |
        Q(can_edit_groups__users__user=user) |
        Q(can_edit_users__user=user)
    ).distinct()


def permitted_to_delete(queryset, user):

    return queryset.filter(Q(owner=user)).distinct()


