# Generated by Django 5.0.2 on 2024-02-09 03:18

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("published", "0010_measurement_access_measurement_owner_and_more"),
    ]

    operations = [
        migrations.RenameField(
            model_name="measurement",
            old_name="sigfig",
            new_name="precision",
        ),
    ]
