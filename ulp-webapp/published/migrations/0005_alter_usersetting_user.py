# Generated by Django 5.1.6 on 2025-04-05 02:49

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("published", "0004_usersetting"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.AlterField(
            model_name="usersetting",
            name="user",
            field=models.OneToOneField(
                on_delete=django.db.models.deletion.DO_NOTHING,
                related_name="setting",
                to=settings.AUTH_USER_MODEL,
            ),
        ),
    ]
