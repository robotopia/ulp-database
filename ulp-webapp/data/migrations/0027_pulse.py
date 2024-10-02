# Generated by Django 5.0.8 on 2024-09-30 12:12

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("data", "0026_workingephemeris_dm"),
    ]

    operations = [
        migrations.CreateModel(
            name="Pulse",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                (
                    "mjd_start",
                    models.FloatField(
                        help_text="The (topocentric) MJD defining the start of the pulse"
                    ),
                ),
                (
                    "mjd_end",
                    models.FloatField(
                        help_text="The (topocentric) MJD defining the end of the pulse"
                    ),
                ),
                (
                    "tags",
                    models.CharField(
                        blank=True,
                        help_text='Comma-separated tags that can be used to categorise different pulses into groups, e.g. "MP", "IP".',
                        max_length=255,
                        null=True,
                    ),
                ),
                (
                    "lightcurve",
                    models.ForeignKey(
                        help_text="The lightcurve to which this pulse belongs.",
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="pulses",
                        to="data.lightcurve",
                    ),
                ),
            ],
            options={
                "ordering": ["lightcurve", "mjd_start"],
            },
        ),
    ]
