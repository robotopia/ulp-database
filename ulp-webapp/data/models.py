from django.db import models
from django.contrib.auth.models import User, Group
from django.core.exceptions import ValidationError
from django.urls import reverse
from published import models as published_models
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation
from astropy.time import Time
import numpy as np
from scipy.stats import vonmises, exponnorm
from scipy.special import erf, erfc, erfcx
from scipy.optimize import curve_fit
from scipy.signal import convolve
from decimal import Decimal

from common.models import AbstractPermission
from common.utils import barycentre, scale_to_frequency


class Observation(AbstractPermission):

    obsid = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="Observation ID",
        verbose_name="ObsID",
    )

    telescope_name = models.CharField(
        max_length=255,
        null=True,
        blank=True,
        help_text="The telescope that made this observation. Must match a string in AstroPy's EarthLocation.get_site_names().",
    )

    freq = models.FloatField(
        null=True,
        blank=True,
        help_text="The centre frequency of this observation, in MHz.",
        verbose_name="Frequency (MHz)",
    )

    bw = models.FloatField(
        null=True,
        blank=True,
        help_text="The bandwidth of this observation, in MHz.",
        verbose_name="Bandwidth (MHz)",
    )

    start_mjd = models.DecimalField(
        decimal_places=20,
        max_digits=30,
        help_text="The MJD of the start of the observation.",
        verbose_name="Start MJD",
    )

    duration = models.FloatField(
        null=True,
        blank=True,
        help_text="The duration of the observation, in seconds.",
        verbose_name="Duration (s)",
    )

    ulps = models.ManyToManyField(
        published_models.Ulp,
        blank=True,
        help_text="ULPs which might be detectable in this observation.",
        related_name="observations",
        verbose_name="ULPs",
    )

    pointing_raj_deg = models.FloatField(
        null=True,
        blank=True,
        help_text="The right ascension (J2000) of the pointing centre, in deg.",
        verbose_name="RA (J2000) (°)",
    )

    pointing_decj_deg = models.FloatField(
        null=True,
        blank=True,
        help_text="The declination (J2000) of the pointing centre, in deg.",
        verbose_name="Dec (J2000) (°)",
    )

    @property
    def start_gps(self):
        t = Time(self.start_mjd, scale='utc', format='mjd')
        return float(t.gps)

    def __str__(self):
        return f'{self.telescope_name} ({self.start_mjd})'

    def clean(self):
        # Make sure the telescope is listed in Astropy
        if self.telescope_name and self.telescope_name not in EarthLocation.get_site_names():
            raise ValidationError(f"\"{self.telescope_name}\" is not found among AstroPy's EarthLocation's site names.")

        # The following checks only apply to new observations
        if self.pk is None:
            # Check for unique Observation IDs (per telescope)
            # Have to do this here (instead of as a unique constraint) because we're allowing
            # these fields to be NULL, in which case the constraint doesn't apply
            if self.telescope_name and self.obsid:
                qs = Observations.objects.filter(telescope_name=self.telescope_name, obsid=self.obsid)
                if qs.exists():
                    raise ValidationError(f"An observation for {self.telescope} already exists with the ObsID {self.obsid}")

            # Check for observations with the same telescope that overlap in time with this one
            if self.telescope_name and self.duration:
                qs = Observations.objects.filter(
                         telescope=self.telescope,
                         start_mjd__gte=this.start_mjd,
                         start_mjd__lte=this.start_mjd + Decimal(this.duration/86400.0),
                     )
                if qs.exists():
                    raise ValidationError(f"At least once observation already exists at this time ({qs.first().start_mjd})")

    class Meta:
        ordering = ['start_mjd', 'freq', 'telescope_name']


class Nondetection(models.Model):

    observation = models.ForeignKey(
        "Observation",
        on_delete=models.CASCADE,
        help_text="The observation being searched for a detection.",
        related_name="nondetections",
    )

    ulp = models.ForeignKey(
        published_models.Ulp,
        on_delete=models.CASCADE,
        help_text="The ULP being searched for.",
        verbose_name="ULP",
        related_name="nondetections",
    )

    local_rms = models.FloatField(
        null=True,
        blank=True,
        help_text="The local rms of the observation, in Jy, which may be used to define an upper limit for the non-detection.",
    )

    def __str__(self):
        return f"Non-detection of {self.ulp} in {self.observation}"

    class Meta:
        verbose_name = "Non-detection"
        verbose_name_plural = "Non-detections"
        constraints = [
            models.UniqueConstraint(fields=['observation', 'ulp'], name="unique_obs_ulp_for_nondetection"),
        ]


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
        verbose_name="Bandwidth (MHz)",
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

    fluence_Jy_s = models.FloatField(
        null=True,
        blank=True,
        help_text="The fluence of the pulse in Jy s.",
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

    def dispersion_offset(self, dm):
        '''
        This calculates the time offset needed to get from the recorded times to the times
        that would have been recorded if the supplied DM was used instead:

          timestamps_in_database + return_value_of_this_function = new_timestamps_assuming_new_dm

        To get the offset needed to recover the timestamps as originally recorded at the telescope,
        call this function with the value of dm = 0.

        The returned value is in units of days, since the times are in MJD.
        '''
        if dm is None:
            return 0.0

        D = 4.148808e3/86400 # Dispersion constant in the appropriate units
        orig_dm_freq = self.dm_freq or np.inf # Recall an empty field (i.e. null value) "means" infinite frequency

        # Offset to get back to to original timestamps
        orig_offset = D * self.dm * (1/self.freq**2 - 1/orig_dm_freq**2)

        # Offset to get to new DM
        new_offset = D * dm / self.freq**2

        # Total offset
        total_offset = orig_offset - new_offset

        return total_offset

    def times(self, pol=None, dm=None):
        if pol is None:
            sample_numbers = np.array([p.sample_number for p in self.points.all()])
        else:
            sample_numbers = np.array([p.sample_number for p in self.points.all() if p.pol == pol])

        times = self.t0 + (sample_numbers * self.dt / 86400)  # Cheaper than astropy units

        # Get the times appropriate for the given DM
        if dm is not None:
            times += self.dispersion_offset(dm=dm)

        return times

    def bary_times(self, pol=None, dm=None):
        times = barycentre(
            self.ulp,
            self.times(pol=pol, dm=dm),
            EarthLocation.of_site(self.telescope),
        )
        return times

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

    ra = models.CharField(
        max_length=31,
        null=True,
        blank=True,
        default="00:00:00",
        help_text="The right ascension (HH:MM:SS.S)",
        verbose_name="RA",
    )

    dec = models.CharField(
        max_length=31,
        null=True,
        blank=True,
        default="00:00:00",
        help_text="The declination (±DD:MM:SS.S)",
    )

    pepoch = models.FloatField(
        null=True,
        blank=True,
        default=60000.0,
        help_text="The reference epoch for period measurement (MJD)",
        verbose_name="PEPOCH",
    )

    p0 = models.FloatField(
        null=True,
        blank=True,
        default=3600.0,
        help_text="The period at time \"pepoch\" (s)",
        verbose_name="Period",
    )

    p1 = models.FloatField(
        null=True,
        blank=True,
        help_text="The period derivative at time \"pepoch\" (s/s)",
        verbose_name="Period derivative",
    )

    pb = models.FloatField(
        null=True,
        blank=True,
        help_text="The orbital period (h)",
        verbose_name="Orbital period",
    )

    dm = models.FloatField(
        null=True,
        blank=True,
        default=0.0,
        help_text="The dispersion measure (pc/cm^3)",
        verbose_name="DM",
    )

    spec_alpha = models.FloatField(
        null=True,
        blank=True,
        help_text="The spectral index. The 'alpha' in the equation S = S1GHz * f^alpha * exp(q*ln(f)^2), where f is the frequency in GHz and S is in Jy.",
        verbose_name="Spectral index",
    )

    spec_q = models.FloatField(
        null=True,
        blank=True,
        help_text="The 'q' in the equation S = S1GHz * f^alpha * exp(q*ln(f)^2), where f is the frequency in GHz and S is in Jy",
        verbose_name="Spectral curvature",
    )

    tausc_1GHz = models.FloatField(
        null=True,
        blank=True,
        help_text="The scattering timescale at 1GHz, in seconds",
        verbose_name="tau_sc",
    )

    covariance = models.OneToOneField(
        "WorkingEphemerisCovariance",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        help_text="The covariance matrix for the ephemeris parameters",
    )

    @property
    def coord(self):
        if self.ra is None or self.dec is None:
            return None

        try:
            coord = SkyCoord(f"{self.ra} {self.dec}", unit=(u.hourangle, u.degree), frame='icrs')
        except:
            return None

        return coord

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

    def fold(self, mjds, freqs_MHz=None):
        # WARNING: Duplicated code!
        # Compare fold() in common.js
        pepoch = self.pepoch
        period = self.p0

        # If freqs_MHz is given, dedisperse to infinite frequency
        if freqs_MHz is not None:
            D = 4.148808e3/86400 # Dispersion constant in the appropriate units
            mjds -= D * self.dm / freqs_MHz**2

        pulses_phases = (mjds - pepoch) / (period/86400.0)
        pulses, phases = np.divmod(pulses_phases + 0.5, 1)
        phases -= 0.5

        return pulses, phases

    def unfold(self, pulses, phases, freqs_MHz=None):
        pepoch = self.pepoch
        period = self.p0

        # If freqs_MHz is given, undedisperse from infinite frequency to given frequencies
        if freqs_MHz is not None:
            D = 4.148808e3/86400 # Dispersion constant in the appropriate units
            mjds += D * self.dm / freqs_MHz**2

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

    def residuals_for_toa_qs(self, toa_qs):
        '''
        Returns the residual in units of days for the ToAs in the supplied queryset
        '''
        mjds = np.array([float(toa.toa_mjd) for toa in toa_qs])
        freqs_MHz = np.array([toa.pulse.lightcurve.freq for toa in toa_qs])
        _, phase = self.fold(mjds, freqs_MHz)
        return phase*self.p0/86400.0

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
        if len(string) > 0:
            array = np.fromstring(string, sep=',')
            return WorkingEphemerisCovariance.from_array(array)
        return

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
        help_text="The (topocentric, not-dedispersed) MJD defining the start of the pulse",
    )

    mjd_end = models.FloatField(
        help_text="The (topocentric,not-dedispersed) MJD defining the end of the pulse",
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

    name = models.CharField(
        max_length=127,
        null=True,
        blank=True,
        help_text="An optional name to help the user distinguish multiple templates for the same source.",
    )

    def values(self, times, tausc=None, freq=None, bw=None, sc_idx=-4.0, nchan=100):
        '''
        Returns template values
        TODO: figure out a sensible normalisation for arbitrary templates
        '''
        sum_of_components = np.sum([component.values(times) for component in self.components.all()], axis=0)

        # Make it a 2D array so that it can be convolved

        # TODO: Normalisation
        #sum_of_weights = np.sum([component.weight for component in self.components.all()]) # <-- bad normalisation?

        if tausc is not None and tausc > 0.0:
            # If freq and bw are not None, then create nchan kernels
            if freq is not None and bw is not None:
                df = bw/nchan
                flo = freq - bw/2
                freqs = np.arange(nchan)*df + flo + df/2
                tauscs = tausc * (freq/freqs).decompose()**sc_idx
            else:
                tauscs = np.array([tausc])
            τ, Times = np.meshgrid(tauscs, times*86400.0)

            # Scattered pulses, when the pulses are gaussians, are equal to exponentially modified Gaussians:
            #   https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
            dynspec = np.zeros(τ.shape)
            dt = (times[1] - times[0])*86400.0
            for component in self.components.all():
                A, μ, σ = component.weight, component.mu*86400.0, component.sigma*86400.0
                dynspec += A*exponnorm.pdf(Times, τ/σ, loc=μ, scale=σ) * np.sqrt(2*np.pi)*σ
            lightcurve = np.mean(dynspec, axis=1)
        else:
            lightcurve = np.squeeze(sum_of_components)

        return lightcurve

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
        res = f"Template ({self.pk}) for {self.working_ephemeris.ulp}"
        if self.name is not None:
            res += f" ({self.name})"
        return res

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
    mu = models.FloatField(help_text="In units of days.")
    sigma = models.FloatField(help_text="In units of days (must be > 0).")

    def values(self, times):
        # If using Gaussians:
        return self.weight*np.exp(-0.5*((times - self.mu)/self.sigma)**2)

        # von Mises function, with everything multiplied by 2pi to convert
        # it to the support expected by the pdf:
        #return self.weight*vonmises.pdf(2*np.pi*phases, loc=2*np.pi*self.mu, kappa=1/(2*np.pi*self.sigma)**2)

    def clean(self):
        if self.sigma <= 0.0:
            raise ValidationError(f"The component cannot have negative width")



class Toa(models.Model):

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    pulse = models.ForeignKey(
        "Pulse",
        on_delete=models.RESTRICT,
        null=True,
        blank=True,
        help_text="The pulse being fitted to.",
        related_name="toas",
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
        help_text="The barycentered (and not dedispersed) MJD of the ToA.",
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
        help_text="The baseline level fitted to the lightcurve (in units of Jy).",
    )

    baseline_slope = models.FloatField(
        null=True,
        blank=True,
        help_text="The baseline slope fitted to the lightcurve (in units of Jy/day).",
    )

    include_in_fit = models.BooleanField(
        default=True,
        help_text="Whether to include this ToA in the timing fit.",
    )

    def get_toa_range(self):
        # Sample the template at N points between the first and last data point
        # "First" point is either the first point in the lightcurve,
        #   or the point closest to -0.5 phase, whichever is latest.
        # "Last" point is either the last point in the lightcurve,
        #   or the point closest to +0.5 phase, whichever is earliest.
        lc = self.pulse.lightcurve
        times = lc.bary_times(dm=0.0)

        half_period = self.template.working_ephemeris.p0/2
        toa_mjd = self.toa_mjd or 0.5*np.sum(self.pulse.bary_start_end())
        pulse_start = float(toa_mjd) - half_period/86400
        pulse_end = float(toa_mjd) + half_period/86400

        lc_start = np.min(times)
        lc_end = np.max(times)

        min_t = pulse_start if lc_start < pulse_start else lc_start
        max_t = pulse_end if lc_end > pulse_end else lc_end

        return min_t, max_t

    @property
    def residual(self):
        '''
        Returns the residual in units of phase
        '''
        we = self.template.working_ephemeris
        _, phase = we.fold(float(self.toa_mjd), self.pulse.lightcurve.freq)
        return phase

    @property
    def next_view(self):
        next_toa = Toa.objects.filter(toa_mjd__gt=self.toa_mjd, template__working_ephemeris=self.template.working_ephemeris).earliest('toa_mjd')
        return reverse('toa_view', args=[next_toa.pk]) if next_toa else ''

    @property
    def prev_view(self):
        prev_toa = Toa.objects.filter(toa_mjd__lt=self.toa_mjd, template__working_ephemeris=self.template.working_ephemeris).latest('toa_mjd')
        return reverse('toa_view', args=[prev_toa.pk]) if prev_toa else ''

    @property
    def out_of_date(self):
        return self.updated < self.template.updated

    @property
    def baseline_degree(self):
        if self.baseline_slope is not None:
            return 1

        if self.baseline_level is not None:
            return 0

        # Otherwise, return None, which means no baseline of any kind is fitted
        return None

    def refit(self, save=True, tausc=None, sc_idx=-4.0, **kwargs):
        '''
        Because of the awkwardness of dealing with lightcurves with generally different
        sampling rates, we simply fit the template to the points, rather than
        trying to use a method based on cross-correlation.

        Keep in mind that the template defines a *shape*, so any fitting function should have
        an "amplitude" as a free parameter. This fitted amplitude can be stored with the ToA
        as a fitted parameter.

        curve_fit doesn't like it when MJDs are used, since the times to be fitted are
        very low down in precision. Better to do everything in terms of time since the start
        and then add it all back afterwards. That's why t0 is preserved (below).

        The default behaviour is to use te existing ToA parameters as the initial guess to
        the fit (i.e. the 'p0' values provided to curve_fit). However, these can be overridden
        by supplying keyword arguments (kwargs) matching the property names, with the values
        to be overridden, e.g., toa.refit(toa_mjd=60000).

        The keywords all take numerical values. In addition, the keywords 'toa_mjd' and 'ampl'
        may also be set to the special value of 'peak', in which case it is set to the MJD
        corresponding to the largest value in the lightcurve. If baseline_level or
        baseline_slope are explicitly set to None, then those parameters are not included
        in the fit.

        tausc is a scattering timescale in seconds. If it is provided, the values are fit to
        the scattered template instead of the "raw" template itself.

        This function does not return anything: it is useful for the side-effect of changing
        the object's fields' values to the fitted values. If the save parameter is True
        (default), the object will be committed to the database.
        '''

        # Get lightcurve for this pulse
        lc = self.pulse.lightcurve
        times = lc.bary_times(dm=0.0)
        values = lc.values()

        min_t, max_t = self.get_toa_range()

        on_pulse = np.logical_and(times >= min_t, times <= max_t)

        times = times[on_pulse]
        values = values[on_pulse]
        t0 = times[0]

        # Set up the initial values
        max_idx = np.nanargmax(values) # For possible use if either 'toa_mjd' or 'ampl' is initialised to 'peak'
        if 'toa_mjd' in kwargs:
            self.toa_mjd = times[max_idx] if kwargs['toa_mjd'] == 'peak' else kwargs['toa_mjd']

        if 'ampl' in kwargs:
            self.ampl = values[max_idx] if kwargs['ampl'] == 'peak' else kwargs['ampl']

        if 'baseline_level' in kwargs:
            self.baseline_level = kwargs['baseline_level']

        if 'baseline_slope' in kwargs:
            self.baseline_slope = kwargs['baseline_slope']

        # Define the most general function to be fitted...
        def template_func(time, toa_mjd, ampl, baseline_level, baseline_slope):
            return ampl*self.template.values(time - toa_mjd,
                    tausc=tausc,
                    freq=lc.freq,
                    bw=lc.bw,
                    sc_idx=sc_idx,) + baseline_level + baseline_slope*(time - toa_mjd)

        # ... but choose the version of this according to which parameters are actually being fitted
        if self.baseline_level is not None and self.baseline_slope is not None:
            p0 = (float(self.toa_mjd) - t0, self.ampl, self.baseline_level, self.baseline_slope)
            bounds = ((-np.inf, 0.0, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf))
            fit_func = template_func
        elif self.baseline_level is not None: # and baseline_slope *is* None
            p0 = (float(self.toa_mjd) - t0, self.ampl, self.baseline_level)
            bounds = ((-np.inf, 0.0, -np.inf), (np.inf, np.inf, np.inf))
            def fit_func(time, toa_mjd, ampl, baseline_level):
                return template_func(time, toa_mjd, ampl, baseline_level, 0.0)
        elif self.baseline_slope is not None: # and baseline_level *is* None
            p0 = (float(self.toa_mjd) - t0, self.ampl, self.baseline_slope)
            bounds = ((-np.inf, 0.0, -np.inf), (np.inf, np.inf, np.inf))
            def fit_func(time, toa_mjd, ampl, baseline_slope):
                return template_func(time, toa_mjd, ampl, 0.0, baseline_slope)
        else: # both baseline_level and baseline_slope are None
            p0 = (float(self.toa_mjd) - t0, self.ampl)
            bounds = ((-np.inf, 0.0), (np.inf, np.inf))
            def fit_func(time, toa_mjd, ampl):
                return template_func(time, toa_mjd, ampl, 0.0, 0.0)

        # Do the fit itself
        popt, pcov = curve_fit(fit_func, times - t0, values, p0=p0, bounds=bounds)

        # Unpack the fitted values. At the moment, we are throwing away the covariances of
        # the two baseline parameters.
        if self.baseline_level is not None and self.baseline_slope is not None:
            self.toa_mjd, self.ampl, self.baseline_level, self.baseline_slope = popt
        elif self.baseline_level is not None: # and baseline_slope *is* None
            self.toa_mjd, self.ampl, self.baseline_level = popt
            self.baseline_slope = None
        elif self.baseline_slope is not None: # and baseline_level *is* None
            self.toa_mjd, self.ampl, self.baseline_slope = popt
            self.baseline_level = None
        else: # both baseline_level and baseline_slope are None
            self.toa_mjd, self.ampl = popt
            self.baseline_level = None
            self.baseline_slope = None

        # Restore the previously subtracted t0 value from the fitted MJD
        self.toa_mjd += t0

        # Unpack the covariances. At the moment, we are throwing away the covariances of
        # the two baseline parameters.
        toa_mjd_err, self.ampl_err = np.sqrt(np.diag(pcov)[:2])
        self.toa_err_s = toa_mjd_err * 86400.0

        if self.ampl_ref_freq is None:
            self.ampl_ref_freq = 1000.0 # This field is (or should be deprecated, but is required in the meantime

        if save:
            self.save()

    def __str__(self):
        return f"ToA from {self.template} ({self.toa_mjd})"

    class Meta:
        verbose_name = "ToA"
        verbose_name_plural = "ToAs"
        ordering = ["toa_mjd"]
