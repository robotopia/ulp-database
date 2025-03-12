# Generated by Django 5.1.6 on 2025-03-12 08:15

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0057_observation_obsid_observation_pointing_decj_deg_and_more'),
        ('published', '0003_ulp_whitelist_users'),
    ]

    operations = [
        migrations.CreateModel(
            name='Nondetection',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('local_rms', models.FloatField(blank=True, help_text='The local rms of the observation, in Jy, which may be used to define an upper limit for the non-detection.', null=True)),
                ('observation', models.ForeignKey(help_text='The observation being searched for a detection.', on_delete=django.db.models.deletion.CASCADE, related_name='nondetections', to='data.observation')),
                ('ulp', models.ForeignKey(help_text='The ULP being searched for.', on_delete=django.db.models.deletion.CASCADE, related_name='nondetections', to='published.ulp', verbose_name='ULP')),
            ],
            options={
                'verbose_name': 'Non-detection',
                'verbose_name_plural': 'Non-detections',
                'constraints': [models.UniqueConstraint(fields=('observation', 'ulp'), name='unique_obs_ulp_for_nondetection')],
            },
        ),
    ]
