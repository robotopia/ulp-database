from django.db import models
from django.contrib.auth.models import User, Group
from django.core.exceptions import ValidationError
from django.urls import reverse
from published import models as published_models
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation
from astropy.time import Time
import numpy as np
from decimal import Decimal

from common.models import AbstractPermission
from common.utils import barycentre

class TimeOfArrival(AbstractPermission):

    ulp = models.ForeignKey(
        published_models.Ulp,
        on_delete=models.CASCADE,
        help_text="The associated ULP.",
        verbose_name="ULP",
        related_name="times_of_arrival",
    )

    mjd = models.DecimalField(
        decimal_places=20,
        max_digits=30,
        help_text="The barycentered MJD of the time of arrival.",
        verbose_name="MJD",
    )

    mjd_err = models.DecimalField(
        decimal_places=20,
        max_digits=30,
        null=True,
        blank=True,
        help_text="The 1σ uncertainty of the time of arrival (in days).",
        verbose_name="MJD error",
    )

    raw_mjd = models.DecimalField(
        decimal_places=20,
        max_digits=30,
        blank=True,
        null=True,
        help_text="The MJD of the time of arrival, before barycentring or dedispersing.",
        verbose_name="Raw MJD",
    )

    telescope_name = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="The telescope that made this detection. Must match a string in AstroPy's EarthLocation.get_site_names().",
    )

    freq = models.FloatField(
        null=True,
        blank=True,
        help_text="The centre frequency of this detection, in MHz.",
        verbose_name="Frequency (MHz)",
    )

    bw = models.FloatField(
        null=True,
        blank=True,
        help_text="The bandwidth of this detection.",
        verbose_name="Bandwidth",
    )

    spectral_index = models.FloatField(
        null=True,
        blank=True,
        help_text="The spectral index of this detection.",
    )

    rotation_measure = models.FloatField(
        null=True,
        blank=True,
        help_text="The rotation measure of this detection in rad/m^2.",
    )

    peak_flux_Jy = models.FloatField(
        null=True,
        blank=True,
        help_text="The peak flux of the pulse in Jy.",
    )

    upper_limit = models.BooleanField(
        default=False,
        help_text="Whether the given peak flux is an upper limit (thereby signifying that this ToA is a \"non-detection\").",
    )

    pulse_width = models.FloatField(
        null=True,
        blank=True,
        help_text="The width of the pulse in seconds.",
    )

    barycentred = models.BooleanField(
        default=True,
        help_text="Whether this TOA has been barycentred.",
    )

    dedispersed = models.BooleanField(
        default=True,
        help_text="Whether this TOA has been corrected for dedispersion.",
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text="Any extra information about the ToA, its measurement method, or the pulse in general.",
    )

    plots = models.ManyToManyField(
        "Plot",
        blank=True,
        help_text="Plots relevant to this ToA.",
        related_name="times_of_arrival",
    )

    lightcurve = models.ForeignKey(
        "Lightcurve",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        help_text="The lightcurve from which this ToA was derived.",
        related_name="toas",
    )

    pulse = models.ForeignKey(
        "Pulse",
        on_delete=models.SET_NULL, # TODO: Change to not-nullable and CASCADE when happy with everything else
        null=True,
        blank=True,
        help_text="The pulse from which this ToA was derived.",
        related_name="toas",
    )

    def __str__(self):
        return f'{self.mjd} ({self.ulp})'

    def clean(self):
        if self.telescope_name and self.telescope_name not in EarthLocation.get_site_names():
            raise ValidationError(f"\"{self.telescope_name}\" is not found among AstroPy's EarthLocation's site names.")

    class Meta:
        ordering = ['ulp', 'mjd']
        verbose_name_plural = "Times of arrival"


class EphemerisParameter(models.Model):

    tempo_name = models.CharField(
        unique=True,
        max_length=63,
        help_text="The name of this parameter according to TEMPO.",
    )

    parameter = models.OneToOneField(
        published_models.Parameter,
        on_delete=models.CASCADE,
        help_text="The parameter type.",
        related_name="ephemeris_parameter",
    )

    tempo_astropy_units = models.CharField(
        max_length=31,
        help_text="The (astropy-conversant) unit if this parameter that TEMPO uses. The unit must be dimensionally equivalent to that of the matching parameter.",
    )

    def validate_unique(self, *args, **kwargs):
        '''
        Ensure that the given astropy units are congruent to the parameter's units.
        '''
        try:
            unit = u.Unit(self.tempo_astropy_units)
        except:
            raise ValidationError(f'{self.tempo_astropy_units} is not a valid AstroPy unit')

        if not unit.is_equivalent(self.parameter.astropy_unit):
            raise ValidationError(f"The units must be equivalent to that of {self.parameter}")

    def save(self, *args, **kwargs):
        self.validate_unique()
        super().save(*args, **kwargs)

    def __str__(self):
        return self.tempo_name

    class Meta:
        ordering = ["tempo_name",]


class EphemerisMeasurement(models.Model):

    ephemeris_parameter = models.ForeignKey(
        "EphemerisParameter",
        on_delete=models.CASCADE,
        help_text="The associated parameter of this measurement.",
        related_name = "ephemeris_measurements",
    )

    measurement = models.ForeignKey(
        published_models.Measurement,
        on_delete=models.CASCADE,
        help_text="The measurement itself.",
        related_name = "ephemeris_measurements",
    )

    @property
    def value(self):
        return self.measurement.astropy_quantity.to(self.ephemeris_parameter.tempo_astropy_units).value

    @property
    def ulp(self):
        return self.measurement.ulp

    @property
    def owner(self):
        return self.measurement.owner

    def validate_unique(self, *args, **kwargs):
        '''
        Custom validation is needed because the unique constraints are too complex.

        First, we must ensure that the measurement is of the same parameter type as
        the ephemeris_parameter.

        Second, we must ensure that each user (measurement.owner) only gets
        one ephemeris per ulp.
        '''
        if self.ephemeris_parameter.parameter != self.measurement.parameter:
            print(f'{self.ephemeris_parameter} vs {self.measurement.parameter}')
            raise ValidationError(f'Parameter type mismatch: the measurement must be an instance of "{self.ephemeris_parameter.parameter}"')

        # Get the user who owns the measurement, and the associated ulp
        user = self.measurement.owner
        ulp = self.measurement.ulp

        # See if this user already owns an ephemeris measurement of this
        # parameter type for the same ULP
        qs = EphemerisMeasurement.objects.filter(
            measurement__owner=user,
            measurement__ulp=ulp,
            ephemeris_parameter=self.ephemeris_parameter,
        ).exclude(
            pk=self.pk
        )

        if qs.exists():
            raise ValidationError(f"{user}'s {ulp} ephemeris already contains the parameter {self.ephemeris_parameter}.")

    def save(self, *args, **kwargs):
        self.validate_unique()
        super().save(*args, **kwargs)

    def __str__(self):
        return f'{self.ephemeris_parameter}: {self.measurement.astropy_quantity.to(self.ephemeris_parameter.tempo_astropy_units)}'

    class Meta:
        ordering = ['measurement__ulp', 'ephemeris_parameter',]


class Plot(models.Model):

    image = models.ImageField(
        upload_to="plots",
        help_text="The plot itself.",
    )

    description = models.TextField(
        null=True,
        blank=True,
    )

    owner = models.ForeignKey(
        User,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="plots",
    )

    def __str__(self):
        return self.image.name



class Lightcurve(AbstractPermission):

    ulp = models.ForeignKey(
        published_models.Ulp,
        on_delete=models.CASCADE,
        help_text="The source being observed.",
        related_name="lightcurves",
    )

    freq = models.FloatField(
        help_text="The centre frequency of this lightcurve, in MHz.",
        verbose_name="Frequency (MHz)",
    )

    bw = models.FloatField(
        help_text="The bandwidth of this lightcurve, in MHz.",
        verbose_name="Bandwidth (MHz)",
    )

    t0 = models.FloatField(
        help_text="Time of first sample (MJD)",
        verbose_name="Start time (MJD)",
    )

    dt = models.FloatField(
        help_text="Duration of each time bin (s)",
        verbose_name="Time bin duration (s)",
    )

    dm = models.FloatField(
        help_text="The dispersion measure used to make this lightcurve, in pc/cm^3.",
        default=0.0,
        verbose_name="DM (pc/cm^3)",
    )

    dm_freq = models.FloatField(
        null=True,
        blank=True,
        help_text="The reference frequency used when dedispersing, in MHz. Leaving this field blank has a special meaning: it is equivalent to setting the reference frequency to (+ve) infinity.",
    )

    telescope = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="The telescope that made this observation. Must match a string in AstroPy's EarthLocation.get_site_names().",
    )

    def times(self, pol=None):
        if pol is None:
            sample_numbers = np.array([p.sample_number for p in self.points.all()])
        else:
            sample_numbers = np.array([p.sample_number for p in self.points.all() if p.pol == pol])
        return self.t0 + (sample_numbers * self.dt / 86400)  # Cheaper than astropy units

    def bary_times(self, pol=None):
        return barycentre(self.ulp, self.times(pol=pol), EarthLocation.of_site(self.telescope))

    def values(self, pol=None):
        if pol is None:
            return np.array([p.value for p in self.points.all()])
        else:
            return np.array([p.value for p in self.points.all() if p.pol == pol])

    @property
    def next_view(self):
        next_lightcurve = self.ulp.lightcurves.filter(t0__gt=self.t0).earliest('t0')
        return reverse('lightcurve_view', args=[next_lightcurve.pk]) if next_lightcurve else ''

    @property
    def prev_view(self):
        prev_lightcurve = self.ulp.lightcurves.filter(t0__lt=self.t0).latest('t0')
        return reverse('lightcurve_view', args=[prev_lightcurve.pk]) if prev_lightcurve else ''

    def __str__(self) -> str:
        return f"Lightcurve ({self.ulp}, {self.t0})"

    class Meta:
        ordering = ['ulp', 't0']


class LightcurvePoint(models.Model):

    lightcurve = models.ForeignKey(
        "Lightcurve",
        on_delete=models.CASCADE,
        help_text="The lightcurve to which this point belongs.",
        related_name="points",
    )

    sample_number = models.IntegerField()

    pol = models.CharField(
        max_length=15,
        default="I",
        help_text="The polarisation of this point.",
    )

    value = models.FloatField(
        help_text="The value of this point, in Jy.",
    )

    err = models.FloatField(
        null=True,
        blank=True,
        help_text="The uncertainty of this point, in Jy.",
    )

    @property
    def t(self):
        return self.lightcurve.t0 + (self.sample_number * self.lightcurve.dt / 86400)  # Cheaper than astropy units

    def __str__(self) -> str:
        return f"LightcurvePoint {self.sample_number} ({self.lightcurve.ulp}, {self.lightcurve.t0})"

    class Meta:
        ordering = ['lightcurve', 'sample_number']
        constraints = [
            models.UniqueConstraint(fields=['lightcurve', 'sample_number'], name="unique_lightcurve_samples"),
        ]


class WorkingEphemeris(AbstractPermission):

    ulp = models.ForeignKey(
        published_models.Ulp,
        on_delete=models.CASCADE,
        help_text="The source this ephemeris applies to",
        related_name="working_ephemerides",
    )

    pepoch = models.FloatField(
        null=True,
        blank=True,
        help_text="The reference epoch for period measurement (MJD)",
    )

    p0 = models.FloatField(
        null=True,
        blank=True,
        help_text="The period at time \"pepoch\" (s)",
    )

    p1 = models.FloatField(
        null=True,
        blank=True,
        help_text="The period derivative at time \"pepoch\" (s/s)",
    )

    pb = models.FloatField(
        null=True,
        blank=True,
        help_text="The orbital period (h)",
    )

    dm = models.FloatField(
        null=True,
        blank=True,
        help_text="The dispersion measure (pc/cm^3)",
    )

    spec_s1GHz = models.FloatField(
        null=True,
        blank=True,
        help_text="The flux density (Jy) at 1 GHz. The 'S1GHz' in the equation S = S1GHz * f^alpha * exp(q*ln(f)^2), where f is the frequency in GHz and S is in Jy.",
    )

    spec_alpha = models.FloatField(
        null=True,
        blank=True,
        help_text="The spectral index. The 'alpha' in the equation S = S1GHz * f^alpha * exp(q*ln(f)^2), where f is the frequency in GHz and S is in Jy.",
    )

    spec_q = models.FloatField(
        null=True,
        blank=True,
        help_text="The 'q' in the equation S = S1GHz * f^alpha * exp(q*ln(f)^2), where f is the frequency in GHz and S is in Jy",
    )

    def predicted_flux_density(self, freq_MHz):
        lnf = np.log(freq_MHz/1e3)
        return np.exp(self.siA*lnf**2 + self.siB*lnf + self.siC)

    def __str__(self) -> str:
        return f"Working ephemeris for {self.ulp} ({self.owner})"

    class Meta:
        verbose_name_plural = "Working ephemerides"
        constraints = [
            models.UniqueConstraint(fields=['owner', 'ulp'], name="working_ephemeris_ulp_owner_unique"),
        ]


class Pulse(models.Model):

    lightcurve = models.ForeignKey(
        "Lightcurve",
        on_delete=models.CASCADE,
        help_text="The lightcurve to which this pulse belongs.",
        related_name="pulses",
    )

    mjd_start = models.FloatField(
        help_text="The (topocentric) MJD defining the start of the pulse",
    )

    mjd_end = models.FloatField(
        help_text="The (topocentric) MJD defining the end of the pulse",
    )

    tags = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="Comma-separated tags that can be used to categorise different pulses into groups, e.g. \"MP\", \"IP\".",
    )

    peak_value_MJD = models.FloatField(
        null=True,
        blank=True,
        help_text="The ToA (topocentric MJD) of the brightest point in the lightcurve within the range defined by this pulse. This field is not intended to be directly editable by the user, but is to be automatically updated whenever the pulse bounds change. It serves to save recalculating this value everytime it's needed (e.g. for a plot).",
    )

    peak_value_Jy = models.FloatField(
        null=True,
        blank=True,
        help_text="The flux density (Jy) of the brightest point in the lightcurve within the range defined by this pulse. This field is not intended to be directly editable by the user, but is to be automatically updated whenever the pulse bounds change. It serves to save recalculating this value everytime it's needed (e.g. for a plot).",
    )

    def __str__(self) -> str:
        return f"Pulse ({self.mjd_start}-{self.mjd_end}) for {self.lightcurve}"

    def clean(self):
        if self.mjd_start >= self.mjd_end:
            raise ValidationError("mjd_start must be less than mjd_end")

    class Meta:
        ordering = ['lightcurve', 'mjd_start']


