# Generated by Django 5.0.2 on 2024-03-21 00:53

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Telescope',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(help_text='The name of the telescope', max_length=255, unique=True)),
                ('abbr', models.CharField(blank=True, help_text='The abbreviation of the name', max_length=15, null=True, unique=True, verbose_name='Abbreviation')),
                ('latitude_deg', models.DecimalField(blank=True, decimal_places=10, help_text='The latitude (in degrees) of the telescope', max_digits=14, null=True, verbose_name='Latitude (°)')),
                ('longitude_deg', models.DecimalField(blank=True, decimal_places=10, help_text='The longitude (in degrees) of the telescope', max_digits=14, null=True, verbose_name='Longitude (°)')),
            ],
            options={
                'ordering': ['name'],
            },
        ),
        migrations.CreateModel(
            name='Backend',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(help_text='The name of the telescope backend', max_length=255, unique=True)),
                ('freq_ctr', models.FloatField(help_text='The centre frequency of the band', verbose_name='Centre frequency')),
                ('bw', models.FloatField(help_text='The bandwidth', verbose_name='Bandwidth')),
                ('freq_units', models.CharField(default='MHz', help_text='An astropy-conversant unit string that applies to the frequencies.', max_length=31)),
                ('telescope', models.ForeignKey(help_text='The telescope to which this backend belongs', on_delete=django.db.models.deletion.CASCADE, related_name='backends', to='data.telescope')),
            ],
            options={
                'ordering': ['telescope', 'name'],
            },
        ),
    ]
