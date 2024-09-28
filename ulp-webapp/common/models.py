from django.db import models
from django.contrib.auth.models import User, Group

class AbstractPermission(models.Model):

    owner = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        null=True,
        blank=True,
        help_text="The owner of this instance.",
        related_name="%(app_label)s_%(class)s_as_owner",
    )

    can_view_groups = models.ManyToManyField(
        Group,
        blank=True,
        help_text="Groups that are granted view privileges.",
        related_name="%(app_label)s_%(class)s_as_viewer",
    )

    can_edit_groups = models.ManyToManyField(
        Group,
        blank=True,
        help_text="Groups that are granted edit privileges.",
        related_name="%(app_label)s_%(class)s_as_editor",
    )

    can_view_users = models.ManyToManyField(
        User,
        blank=True,
        help_text="Users that are granted view privileges.",
        related_name="%(app_label)s_%(class)s_as_viewer",
    )

    can_edit_users = models.ManyToManyField(
        User,
        blank=True,
        help_text="Users that are granted edit privileges.",
        related_name="%(app_label)s_%(class)s_as_editor",
    )

    class Meta:
        abstract = True
