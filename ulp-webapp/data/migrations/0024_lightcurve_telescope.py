# Generated by Django 5.0.8 on 2024-09-29 14:12

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("data", "0023_alter_lightcurvepoint_pol"),
    ]

    operations = [
        migrations.AddField(
            model_name="lightcurve",
            name="telescope",
            field=models.CharField(
                blank=True,
                help_text="The telescope that made this observation. Must match a string in AstroPy's EarthLocation.get_site_names().",
                max_length=255,
                null=True,
            ),
        ),
    ]
