from django.db import models
from django.contrib.auth.models import User, Group
from django.core.exceptions import ValidationError
from django.urls import reverse
from published import models as published_models
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation
from astropy.time import Time
import numpy as np
from scipy.stats import vonmises
from decimal import Decimal

from common.models import AbstractPermission
from common.utils import barycentre, scale_to_frequency, dm_correction

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

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

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

    @property
    def t0_gps(self):
        return int(np.round(Time(self.t0, scale='utc', format='mjd').gps))

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

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

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

    covariance = models.OneToOneField(
        "WorkingEphemerisCovariance",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        help_text="The covariance matrix for the ephemeris parameters",
    )

    @property
    def pepoch_err(self):
        if self.covariance:
            if self.covariance.pepoch_pepoch:
                return np.sqrt(self.covariance.pepoch_pepoch)

    @property
    def p0_err(self):
        if self.covariance:
            if self.covariance.p0_p0:
                return np.sqrt(self.covariance.p0_p0)

    @property
    def p1_err(self):
        if self.covariance:
            if self.covariance.p1_p1:
                return np.sqrt(self.covariance.p1_p1)

    @property
    def pb_err(self):
        if self.covariance:
            if self.covariance.pb_pb:
                return np.sqrt(self.covariance.pb_pb)

    @property
    def dm_err(self):
        if self.covariance:
            if self.covariance.dm_dm:
                return np.sqrt(self.covariance.dm_dm)

    @property
    def spec_alpha_err(self):
        if self.covariance:
            if self.covariance.spec_alpha_spec_alpha:
                return np.sqrt(self.covariance.spec_alpha_spec_alpha)

    @property
    def spec_q_err(self):
        if self.covariance:
            if self.covariance.spec_q_spec_q:
                return np.sqrt(self.covariance.spec_q_spec_q)

    def predicted_flux_density(self, freq_MHz):
        lnf = np.log(freq_MHz/1e3)
        return np.exp(self.siA*lnf**2 + self.siB*lnf + self.siC)

    def fold(self, mjds):
        # WARNING: Duplicated code!
        # Compare fold() in common.js
        pepoch = self.pepoch
        period = self.p0

        pulses_phases = (mjds - pepoch) / (period/86400.0)
        pulses, phases = np.divmod(pulses_phases + 0.5, 1)
        phases -= 0.5

        return pulses, phases

    def unfold(self, pulses, phases):
        pepoch = self.pepoch
        period = self.p0
        mjds = pepoch + (period/86400.0)*(pulses + phases)
        return mjds

    def pulses_from_pulse_number(self, pulse_number):
        pulse_start, pulse_end = self.unfold(pulse_number, np.array([-0.5, 0.5]))
        # ^^^ These still need to be barycentred
        pulses = [] # <-- where we will put the pulses that match
        for pulse in Pulse.objects.all():
            # ^^^ I can't think of a cheaper way of doing this, since Earth locations
            # can differ from pulse to pulse, and since barycentric conversions can't
            # be done at the database level...
            bary_mjd_start, _ = pulse.bary_start_end() 
            if bary_mjd_start >= pulse_start and bary_mjd_start <= pulse_end:
                pulses.append(pulse)
        return pulses

    def extract_lightcurves_from_pulse_number(self, pulse_number, freq_target_MHz):

        # Get the pulses we'll be working with
        pulses = self.pulses_from_pulse_number(pulse_number)

        if len(pulses) == 0:
            raise Exception(f"No pulses with this pulse number ({pulse_number})")

        # Extract the lightcurves, shift them to infinite frequency, and
        # scale them all to the same (arbitrary) frequency
        times = np.concatenate([p.lightcurve.bary_times() + dm_correction(p.lightcurve, self) for p in pulses])
        values = np.concatenate([scale_to_frequency(
            p.lightcurve.freq,
            p.lightcurve.values(),
            freq_target_MHz,
            self.spec_alpha,
            self.spec_q,
        ) for p in pulses])

        return times, values

    def __str__(self) -> str:
        return f"Working ephemeris for {self.ulp}"

    class Meta:
        verbose_name_plural = "Working ephemerides"
        constraints = [
            models.UniqueConstraint(fields=['owner', 'ulp'], name="working_ephemeris_ulp_owner_unique"),
        ]


class WorkingEphemerisCovariance(models.Model):

    pepoch_pepoch = models.FloatField(
        null=True,
        blank=True,
    )

    pepoch_p0 = models.FloatField(
        null=True,
        blank=True,
    )

    pepoch_p1 = models.FloatField(
        null=True,
        blank=True,
    )

    pepoch_pb = models.FloatField(
        null=True,
        blank=True,
    )

    pepoch_dm = models.FloatField(
        null=True,
        blank=True,
    )

    pepoch_spec_alpha = models.FloatField(
        null=True,
        blank=True,
    )

    pepoch_spec_q = models.FloatField(
        null=True,
        blank=True,
    )

    p0_p0 = models.FloatField(
        null=True,
        blank=True,
    )

    p0_p1 = models.FloatField(
        null=True,
        blank=True,
    )

    p0_pb = models.FloatField(
        null=True,
        blank=True,
    )

    p0_dm = models.FloatField(
        null=True,
        blank=True,
    )

    p0_spec_alpha = models.FloatField(
        null=True,
        blank=True,
    )

    p0_spec_q = models.FloatField(
        null=True,
        blank=True,
    )

    p1_p1 = models.FloatField(
        null=True,
        blank=True,
    )

    p1_pb = models.FloatField(
        null=True,
        blank=True,
    )

    p1_dm = models.FloatField(
        null=True,
        blank=True,
    )

    p1_spec_alpha = models.FloatField(
        null=True,
        blank=True,
    )

    p1_spec_q = models.FloatField(
        null=True,
        blank=True,
    )

    pb_pb = models.FloatField(
        null=True,
        blank=True,
    )

    pb_dm = models.FloatField(
        null=True,
        blank=True,
    )

    pb_spec_alpha = models.FloatField(
        null=True,
        blank=True,
    )

    pb_spec_q = models.FloatField(
        null=True,
        blank=True,
    )

    dm_dm = models.FloatField(
        null=True,
        blank=True,
    )

    dm_spec_alpha = models.FloatField(
        null=True,
        blank=True,
    )

    dm_spec_q = models.FloatField(
        null=True,
        blank=True,
    )

    spec_alpha_spec_alpha = models.FloatField(
        null=True,
        blank=True,
    )

    spec_alpha_spec_q = models.FloatField(
        null=True,
        blank=True,
    )

    spec_q_spec_q = models.FloatField(
        null=True,
        blank=True,
    )

    # A helper function to create a WorkingEphemerisCovariance instance from an array
    def from_array(array):
        # 'array' is expected to be in the format output by as_array()
        return WorkingEphemerisCovariance(
            pepoch_pepoch=array[0] if not np.isnan(array[0]) else None,
            pepoch_p0=array[1] if not np.isnan(array[1]) else None,
            pepoch_p1=array[2] if not np.isnan(array[2]) else None,
            pepoch_pb=array[3] if not np.isnan(array[3]) else None,
            pepoch_dm=array[4] if not np.isnan(array[4]) else None,
            pepoch_spec_alpha=array[5] if not np.isnan(array[5]) else None,
            pepoch_spec_q=array[6] if not np.isnan(array[6]) else None,
            p0_p0=array[7] if not np.isnan(array[7]) else None,
            p0_p1=array[8] if not np.isnan(array[8]) else None,
            p0_pb=array[9] if not np.isnan(array[9]) else None,
            p0_dm=array[10] if not np.isnan(array[10]) else None,
            p0_spec_alpha=array[11] if not np.isnan(array[11]) else None,
            p0_spec_q=array[12] if not np.isnan(array[12]) else None,
            p1_p1=array[13] if not np.isnan(array[13]) else None,
            p1_pb=array[14] if not np.isnan(array[14]) else None,
            p1_dm=array[15] if not np.isnan(array[15]) else None,
            p1_spec_alpha=array[16] if not np.isnan(array[16]) else None,
            p1_spec_q=array[17] if not np.isnan(array[17]) else None,
            pb_pb=array[18] if not np.isnan(array[18]) else None,
            pb_dm=array[19] if not np.isnan(array[19]) else None,
            pb_spec_alpha=array[20] if not np.isnan(array[20]) else None,
            pb_spec_q=array[21] if not np.isnan(array[21]) else None,
            dm_dm=array[22] if not np.isnan(array[22]) else None,
            dm_spec_alpha=array[23] if not np.isnan(array[23]) else None,
            dm_spec_q=array[24] if not np.isnan(array[24]) else None,
            spec_alpha_spec_alpha=array[25] if not np.isnan(array[25]) else None,
            spec_alpha_spec_q=array[26] if not np.isnan(array[26]) else None,
            spec_q_spec_q=array[27] if not np.isnan(array[27]) else None,
        )

    def from_str(string):
        array = np.fromstring(string, sep=',')
        return WorkingEphemerisCovariance.from_array(array)

    def as_matrix(self):
        return np.array(
            [
                self.pepoch_pepoch, self.pepoch_p0, self.pepoch_p1, self.pepoch_pb, self.pepoch_dm, self.pepoch_spec_alpha, self.pepoch_spec_q,
                self.pepoch_p0, self.p0_p0, self.p0_p1, self.p0_pb, self.p0_dm, self.p0_spec_alpha, self.p0_spec_q,
                self.pepoch_p1, self.p0_p1, self.p1_p1, self.p1_pb, self.p1_dm, self.p1_spec_alpha, self.p1_spec_q,
                self.pepoch_pb, self.p0_pb, self.p1_pb, self.pb_pb, self.pb_dm, self.pb_spec_alpha, self.pb_spec_q,
                self.pepoch_dm, self.p0_dm, self.p1_dm, self.pb_dm, self.dm_dm, self.dm_spec_alpha, self.dm_spec_q,
                self.pepoch_spec_alpha, self.p0_spec_alpha, self.p1_spec_alpha, self.pb_spec_alpha, self.dm_spec_alpha, self.spec_alpha_spec_alpha, self.spec_alpha_spec_q,
                self.pepoch_spec_q, self.p0_spec_q, self.p1_spec_q, self.pb_spec_q, self.dm_spec_q, self.spec_alpha_spec_q, self.spec_q_spec_q,
            ]
        )

    def as_array(self):
        # If this function is changed, from_array() must also be updated
        return np.array(
            [
                self.pepoch_pepoch, self.pepoch_p0, self.pepoch_p1, self.pepoch_pb, self.pepoch_dm, self.pepoch_spec_alpha, self.pepoch_spec_q, self.p0_p0, self.p0_p1, self.p0_pb, self.p0_dm, self.p0_spec_alpha, self.p0_spec_q, self.p1_p1, self.p1_pb, self.p1_dm, self.p1_spec_alpha, self.p1_spec_q, self.pb_pb, self.pb_dm, self.pb_spec_alpha, self.pb_spec_q, self.dm_dm, self.dm_spec_alpha, self.dm_spec_q, self.spec_alpha_spec_alpha, self.spec_alpha_spec_q, self.spec_q_spec_q,
            ]
        )

    def __str__(self):
        return ",".join([str(x) if x else 'nan' for x in self.as_array()])


class Pulse(models.Model):

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

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

    def bary_start_end(self):
        # Converts the start and end times to barycentre
        return barycentre(self.lightcurve.ulp, [self.mjd_start, self.mjd_end], EarthLocation.of_site(self.lightcurve.telescope))

    def __str__(self) -> str:
        return f"Pulse ({self.mjd_start}-{self.mjd_end}) for {self.lightcurve}"

    def clean(self):
        if self.mjd_start >= self.mjd_end:
            raise ValidationError("mjd_start must be less than mjd_end")

    class Meta:
        ordering = ['lightcurve', 'mjd_start']


class Template(AbstractPermission):

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    included_pulses = models.ManyToManyField(
        "Pulse",
        through="TemplatePulse",
        through_fields=("template", "pulse"),
        related_name="templates",
    )

    working_ephemeris = models.ForeignKey(
        "WorkingEphemeris",
        on_delete=models.CASCADE,
        help_text="The ephemeris used for constructing this template.",
        related_name="templates",
    )

    def values(self, phases):
        '''
        Returns template values normalised so that the whole template has unity area under curve
        '''
        sum_of_components = np.sum([component.values(phases) for component in self.components.all()], axis=0)
        sum_of_weights = np.sum([component.weight for component in self.components.all()])
        return sum_of_components/sum_of_weights

    def values_npoints(self, npoints, ph_ctr=0.0):
        phases = np.linspace(ph_ctr-0.5, ph_ctr+0.5, num=npoints, endpoint=False)
        return phases, self.values(phases)

    def values_dph(self, dph, ph_ctr=0.0):
        # For a given "delta phase" (dph), it won't generally be the case that a whole number
        # of samples will fit into one period. The strategy is to extend slightly beyond
        # the phase range (ph_ctr-0.5, ph_ctr+0.5) as minimally as possible. It is assumed that
        # contributions near the extremes of this range are minimal
        npoints = int(np.ceil(1/dph))
        phase_range = npoints*dph
        phases = np.linspace(ph_ctr - phase_range/2, ph_ctr + phase_range/2, num=npoints, endpoint=False)
        return phases, self.values(npoints)

    @property
    def out_of_date(self):
        return self.updated < self.working_ephemeris.updated

    def __str__(self) -> str:
        return f"Template for {self.working_ephemeris.ulp}"

    class Meta:
        ordering = ['working_ephemeris', 'updated']


class TemplatePulse(models.Model):

    template = models.ForeignKey(
        "Template",
        on_delete=models.CASCADE,
    )

    pulse = models.ForeignKey(
        "Pulse",
        on_delete=models.CASCADE,
    )

    def clean(self):
        if self.pulse.lightcurve.ulp != self.template.working_ephemeris.ulp:
            raise ValidationError(f"The pulse's Ulp ({self.pulse.lightcurve.ulp}) must be the same as the template's Ulp ({self.template.working_ephemeris.ulp}).")


class TemplateComponent(models.Model):

    template = models.ForeignKey(
        "Template",
        on_delete=models.CASCADE,
        help_text="The template to which this component belongs.",
        related_name="components",
    )

    weight = models.FloatField()
    mu = models.FloatField(help_text="In units of pulse phase.")
    sigma = models.FloatField(help_text="In units of pulse phase.")

    def values(self, phases):
        # If using Gaussians, which we're not:
        #return self.weight*np.exp(0.5*((phases - self.mu)/self.sigma)**2)

        # von Mises function, with everything multiplied by 2pi to convert
        # it to the support expected by the pdf:
        return self.weight*vonmises.pdf(2*np.pi*phases, loc=2*np.pi*self.mu, kappa=1/(2*np.pi*self.sigma)**2)

    def clean(self):
        if self.sigma <= 0.0:
            raise ValidationError(f"The component cannot have negative width")



class Toa(models.Model):

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    pulse_number = models.IntegerField(
        help_text="The pulse number as derived from the linked working ephemeris.",
    )

    template = models.ForeignKey(
        "Template",
        on_delete=models.CASCADE,
        help_text="The template used to derive this ToA.",
        related_name="toas",
    )

    toa_mjd = models.DecimalField(
        decimal_places=20,
        max_digits=30,
        help_text="The barycentered MJD of the ToA.",
        verbose_name="ToA (MJD)",
    )

    toa_err_s = models.FloatField(
        help_text="The 1σ uncertainty of the time of arrival (in seconds).",
        verbose_name="ToA error (s)",
    )

    toa_freq_MHz = models.FloatField(
        null=True,
        blank=True,
        help_text="The frequency to be used for this ToA.",
        verbose_name="Frequency (MHz)",
    )

    ampl = models.FloatField(
        help_text="The amplitude of the fitted template that matched the data scaled to 'ampl_ref_freq'.",
    )

    ampl_err = models.FloatField(
        help_text="The uncertainaty of the amplitude.",
    )

    ampl_ref_freq = models.FloatField(
        help_text="The reference frequency used (when scaling the data) to derive the amplitude.",
    )

    baseline_level = models.FloatField(
        null=True,
        blank=True,
        help_text="The baseline level fitted to the lightcurve.",
    )

    baseline_slope = models.FloatField(
        null=True,
        blank=True,
        help_text="The baseline slope fitted to the lightcurve.",
    )

    include_in_fit = models.BooleanField(
        default=True,
        help_text="Whether to include this ToA in the timing fit.",
    )

    @property
    def residual(self):
        '''
        Returns the residual in units of phase
        '''
        we = self.template.working_ephemeris
        _, phase = we.fold(float(self.toa_mjd))
        return phase

    @property
    def next_view(self):
        next_toa = self.template.toas.filter(toa_mjd__gt=self.toa_mjd).earliest('toa_mjd')
        return reverse('toa_view', args=[next_toa.pk]) if next_toa else ''

    @property
    def prev_view(self):
        prev_toa = self.template.toas.filter(toa_mjd__lt=self.toa_mjd).latest('toa_mjd')
        return reverse('toa_view', args=[prev_toa.pk]) if prev_toa else ''

    @property
    def out_of_date(self):
        return self.updated < self.template.updated

    def __str__(self):
        return f"ToA from {self.template} ({self.toa_mjd})"

    class Meta:
        verbose_name = "ToA"
        verbose_name_plural = "ToAs"
        ordering = ["pulse_number"]
        constraints = [
            models.UniqueConstraint(fields=['template', 'pulse_number'], name="unique_toa_per_pulse_and_template"),
        ]
