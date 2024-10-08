# Generated by Django 5.0.8 on 2024-10-05 05:53

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("data", "0040_remove_toa_working_ephemeris_toa_ampl_toa_ampl_err_and_more"),
    ]

    operations = [
        migrations.AlterModelOptions(
            name="toa",
            options={
                "ordering": ["pk", "pulse_number"],
                "verbose_name": "ToA",
                "verbose_name_plural": "ToAs",
            },
        ),
        migrations.AlterField(
            model_name="toa",
            name="toa_err_s",
            field=models.FloatField(
                help_text="The 1σ uncertainty of the time of arrival (in seconds).",
                verbose_name="ToA error (s)",
            ),
        ),
    ]
