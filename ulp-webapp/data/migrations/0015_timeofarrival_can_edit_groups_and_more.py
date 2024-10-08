# Generated by Django 5.0.8 on 2024-09-28 03:12

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("auth", "0012_alter_user_first_name_max_length"),
        ("data", "0014_timeofarrival_telescope_name"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.AddField(
            model_name="timeofarrival",
            name="can_edit_groups",
            field=models.ManyToManyField(
                blank=True,
                help_text="Groups that are granted edit privileges.",
                related_name="%(app_label)s_%(class)s_as_editor",
                to="auth.group",
            ),
        ),
        migrations.AddField(
            model_name="timeofarrival",
            name="can_edit_users",
            field=models.ManyToManyField(
                blank=True,
                help_text="Users that are granted edit privileges.",
                related_name="%(app_label)s_%(class)s_as_editor",
                to=settings.AUTH_USER_MODEL,
            ),
        ),
        migrations.AddField(
            model_name="timeofarrival",
            name="can_view_groups",
            field=models.ManyToManyField(
                blank=True,
                help_text="Groups that are granted view privileges.",
                related_name="%(app_label)s_%(class)s_as_viewer",
                to="auth.group",
            ),
        ),
        migrations.AddField(
            model_name="timeofarrival",
            name="can_view_users",
            field=models.ManyToManyField(
                blank=True,
                help_text="Users that are granted view privileges.",
                related_name="%(app_label)s_%(class)s_as_viewer",
                to=settings.AUTH_USER_MODEL,
            ),
        ),
        migrations.AddField(
            model_name="timeofarrival",
            name="owner",
            field=models.ForeignKey(
                blank=True,
                help_text="The owner of this instance.",
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="%(app_label)s_%(class)s_as_owner",
                to=settings.AUTH_USER_MODEL,
            ),
        ),
    ]
