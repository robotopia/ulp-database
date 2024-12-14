# Generated by Django 5.0.10 on 2024-12-14 04:08

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0051_alter_toa_baseline_level_alter_toa_baseline_slope'),
    ]

    operations = [
        migrations.AddField(
            model_name='workingephemeris',
            name='tausc_1GHz',
            field=models.FloatField(blank=True, help_text='The scattering timescale at 1GHz, in seconds', null=True, verbose_name='τ_{sc,1GHz} (s)'),
        ),
    ]
