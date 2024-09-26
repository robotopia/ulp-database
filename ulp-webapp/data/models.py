from django.db import models
from django.contrib.auth.models import User, Group
from django.core.exceptions import ValidationError
from published import models as published_models
import astropy.units as u
from astropy.coordinates import Angle
from decimal import Decimal

# Create your models here.
class Telescope(models.Model):

    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="The name of the telescope.",
    )

    abbr = models.CharField(
        max_length=15,
        unique=True,
        blank=True,
        null=True,
        help_text="The abbreviation of the name.",
        verbose_name="Abbreviation",
    )

    latitude_deg = models.DecimalField(
        decimal_places=10,
        max_digits=14,
        null=True,
        blank=True,
        help_text="The latitude (in degrees) of the telescope.",
        verbose_name="Latitude (°)",
    )

    longitude_deg = models.DecimalField(
        decimal_places=10,
        max_digits=14,
        null=True,
        blank=True,
        help_text="The longitude (in degrees) of the telescope.",
        verbose_name="Longitude (°)",
    )

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name',]

class Backend(models.Model):

    telescope = models.ForeignKey(
        "Telescope",
        on_delete=models.CASCADE,
        help_text="The telescope to which this backend belongs.",
        related_name="backends",
    )

    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="The name of the telescope backend",
    )

    freq_ctr = models.FloatField(
        null=True,
        blank=True,
        help_text="The centre frequency of the band.",
        verbose_name="Centre frequency",
    )

    bw = models.FloatField(
        null=True,
        blank=True,
        help_text="The bandwidth.",
        verbose_name="Bandwidth",
    )

    freq_units = models.CharField(
        max_length=31,
        null=True,
        blank=True,
        default="MHz",
        help_text="An astropy-conversant unit string that applies to the frequencies.",
    )

    def __str__(self):
        return f'{self.name} ({self.telescope})'

    class Meta:
        ordering = ['telescope', 'name',]


class TimeOfArrival(models.Model):

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

    telescope = models.ForeignKey(
        "Telescope",
        on_delete=models.DO_NOTHING,
        null=True,
        blank=True,
        help_text="The telescope that made this detection.",
        related_name="toas",
    )

    backend = models.ForeignKey(
        "Backend",
        on_delete=models.DO_NOTHING,
        null=True,
        blank=True,
        help_text="The telescope backend that made this detection.",
        related_name="toas",
    )

    freq = models.FloatField(
        null=True,
        blank=True,
        help_text="The centre frequency of this detection.",
        verbose_name="Frequency",
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

    freq_units = models.CharField(
        max_length=31,
        null=True,
        blank=True,
        default="MHz",
        help_text="An astropy-conversant unit string that applies to the frequency.",
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

    def __str__(self):
        return f'{self.mjd} ({self.ulp})'

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


