# Generated by Django 5.0.11 on 2025-03-04 06:50

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0054_workingephemeris_dec_workingephemeris_ra'),
        ('published', '0003_ulp_whitelist_users'),
    ]

    operations = [
        migrations.AlterField(
            model_name='timeofarrival',
            name='bw',
            field=models.FloatField(blank=True, help_text='The bandwidth of this detection.', null=True, verbose_name='Bandwidth (MHz)'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='dec',
            field=models.CharField(blank=True, default='00:00:00', help_text='The declination (±DD:MM:SS.S)', max_length=31, null=True),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='dm',
            field=models.FloatField(blank=True, default=0.0, help_text='The dispersion measure (pc/cm^3)', null=True, verbose_name='DM'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='p0',
            field=models.FloatField(blank=True, default=3600.0, help_text='The period at time "pepoch" (s)', null=True, verbose_name='Period'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='p1',
            field=models.FloatField(blank=True, help_text='The period derivative at time "pepoch" (s/s)', null=True, verbose_name='Period derivative'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='pb',
            field=models.FloatField(blank=True, help_text='The orbital period (h)', null=True, verbose_name='Orbital period'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='pepoch',
            field=models.FloatField(blank=True, default=60000.0, help_text='The reference epoch for period measurement (MJD)', null=True, verbose_name='PEPOCH'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='ra',
            field=models.CharField(blank=True, default='00:00:00', help_text='The right ascension (HH:MM:SS.S)', max_length=31, null=True, verbose_name='RA'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='spec_alpha',
            field=models.FloatField(blank=True, help_text="The spectral index. The 'alpha' in the equation S = S1GHz * f^alpha * exp(q*ln(f)^2), where f is the frequency in GHz and S is in Jy.", null=True, verbose_name='Spectral index'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='spec_q',
            field=models.FloatField(blank=True, help_text="The 'q' in the equation S = S1GHz * f^alpha * exp(q*ln(f)^2), where f is the frequency in GHz and S is in Jy", null=True, verbose_name='Spectral curvature'),
        ),
        migrations.AlterField(
            model_name='workingephemeris',
            name='tausc_1GHz',
            field=models.FloatField(blank=True, help_text='The scattering timescale at 1GHz, in seconds', null=True, verbose_name='tau_sc'),
        ),
        migrations.CreateModel(
            name='Observation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('telescope_name', models.CharField(blank=True, help_text="The telescope that made this observation. Must match a string in AstroPy's EarthLocation.get_site_names().", max_length=255, null=True)),
                ('freq', models.FloatField(blank=True, help_text='The centre frequency of this observation, in MHz.', null=True, verbose_name='Frequency (MHz)')),
                ('bw', models.FloatField(blank=True, help_text='The bandwidth of this observation, in MHz.', null=True, verbose_name='Bandwidth (MHz)')),
                ('start_mjd', models.DecimalField(decimal_places=20, help_text='The MJD of the start of the observation.', max_digits=30, verbose_name='Start MJD')),
                ('duration', models.FloatField(blank=True, help_text='The duration of the observation, in seconds.', null=True, verbose_name='Duration (s)')),
                ('ulps', models.ManyToManyField(blank=True, help_text='ULPs which might be detectable in this observation.', related_name='observations', to='published.ulp', verbose_name='ULPs')),
            ],
            options={
                'ordering': ['start_mjd', 'freq', 'telescope_name'],
            },
        ),
    ]
