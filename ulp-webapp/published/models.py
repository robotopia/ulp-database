from django.db import models
from django.contrib.auth.models import User, Group

from decimal import Decimal

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.constants import c

# Create your models here.
class Ulp(models.Model):

    name = models.CharField(
        max_length=63,
        unique=True,
        help_text="The official name of the ULP",
    )

    abbr = models.CharField(
        max_length=15,
        unique=True,
        blank=True,
        null=True,
        help_text="A short version of the name",
    )

    discoverer = models.ForeignKey(
        Group,
        null=True,
        blank=True,
        help_text="The group who first discovered this ULP",
        on_delete=models.RESTRICT,
        related_name="ulps",
    )

    data_access_groups = models.ManyToManyField(
        Group,
        blank=True,
        help_text="The groups that are allowed access to this ULP's data.",
        related_name="data_accessible_ulps",
    )

    whitelist_users = models.ManyToManyField(
        User,
        blank=True,
        help_text="Any individual users that are given explicit access to this ULP's data.",
        related_name="data_accessible_ulps",
    )

    @property
    def best_progenitor_claims(self):
        claims = self.progenitor_claims.all().order_by('tentative')
        if not claims.exists():
            return None
        least_tentative = claims.first().tentative
        least_tentative_claims = claims.filter(tentative=least_tentative)
        return ', '.join([f'{claim}' for claim in least_tentative_claims])

    class Meta:
        verbose_name = "ULP"
        verbose_name_plural = "ULPs"
        ordering = ['name']

    def __str__(self):
        return self.name


class Article(models.Model):

    citet_text = models.CharField(
        max_length=127,
        blank=True,
        null=True,
        help_text="The display text when this is cited inline, e.g., Smith et al. (1900)",
    )

    title = models.CharField(
        max_length=1023,
        blank=True,
        null=True,
        help_text="The title of the paper",
    )

    doi = models.CharField(
        max_length=1023,
        blank=True,
        null=True,
        help_text="The DOI of the article",
        verbose_name="DOI",
    )

    toa_url = models.CharField(
        max_length=1023,
        blank=True,
        null=True,
        help_text="URL for where ToAs in this paper are published.",
        verbose_name="URL for ToAs",
    )

    eph_url = models.CharField(
        max_length=1023,
        blank=True,
        null=True,
        help_text="URL for where the ephemeris in this paper are published.",
        verbose_name="URL for ephemeris",
    )

    @property
    def doi_url(self):
        return "https://doi.org/" + self.doi

    def __str__(self):
        return self.citet_text

    class Meta:
        ordering = ['citet_text']


class Parameter(models.Model):

    name = models.CharField(
        max_length=63,
        unique=True,
        help_text="The name of the parameter",
    )

    description = models.TextField(
        blank=True,
        null=True,
        help_text="A description of this measurable parameter",
    )

    ascii_symbol = models.CharField(
        max_length=15,
        help_text="The symbol used for this parameter, expressed in ASCII characters",
    )

    latex_symbol = models.CharField(
        max_length=63,
        blank=True,
        null=True,
        help_text="The symbol used for this parameter, expressed as a LaTeX expression",
    )

    astropy_unit = models.CharField(
        max_length=31,
        blank=True,
        null=True,
        help_text="The (astropy-conversant) unit of this parameter (e.g. 's', 'lightyear')",
    )

    def __str__(self):
        if self.unicode_unit is not None:
            return f"{self.name} ({self.unicode_unit})"
        else:
            return f"{self.name}"

    @property
    def unicode_unit(self):
        if self.astropy_unit is None:
            return None

        if u.Unit(self.astropy_unit) == u.dimensionless_unscaled:
            return None

        unicode_unit = u.Unit(self.astropy_unit).to_string(format='unicode')
        if unicode_unit[-1] == '☉':
            unicode_unit = f"{unicode_unit[:-1]}<sub>☉</sub>"

        return unicode_unit

    class Meta:
        ordering = ['name']


class Measurement(models.Model):

    ACCESS_PUBLIC = 'P'
    ACCESS_GROUP = 'G'
    ACCESS_PRIVATE = 'O'

    BODY_TOTAL = 'T'
    BODY_PRIMARY = '1'
    BODY_SECONDARY = '2'
    BODY_TERTIARY = '3'

    BODY_CHOICES=[
        ('T', 'Total'),
        ('1', 'Primary'),
        ('2', 'Secondary'),
        ('3', 'Tertiary'),
    ]

    parameter = models.ForeignKey(
        "Parameter",
        help_text="The parameter beaing measured",
        related_name="measurements",
        on_delete=models.RESTRICT,
    )

    body = models.CharField(
        max_length=2,
        null=True,
        blank=True,
        choices=BODY_CHOICES,
        help_text="If system contains more than one stellar body, which body or bodies this measurement pertains to.",
    )

    ulp = models.ForeignKey(
        "Ulp",
        help_text="The ULP for to which this measurement applies",
        related_name="measurements",
        on_delete=models.RESTRICT,
        verbose_name="ULP",
    )

    quantity = models.DecimalField(
        decimal_places=20,
        max_digits=40,
        help_text="The numerical value of this measurement",
    )

    power_of_10 = models.IntegerField(
        default=0,
        #help_text="The true quantity is multiplied by this power of 10.",
        verbose_name="× 10^",
    )

    precision = models.IntegerField(
        null=True,
        blank=True,
        help_text="The number of decimal places to display",
    )

    err = models.DecimalField(
        decimal_places=20,
        max_digits=40,
        null=True,
        blank=True,
        help_text="The (symmetrical) uncertainty on the quantity.",
        verbose_name="Error (±)",
    )

    err_hi = models.DecimalField(
        decimal_places=20,
        max_digits=40,
        null=True,
        blank=True,
        help_text="The (asymmetrical) upper uncertainty on the quantity.",
        verbose_name="Error high (+)",
    )

    err_lo = models.DecimalField(
        decimal_places=20,
        max_digits=40,
        null=True,
        blank=True,
        help_text="The (asymmetrical) lower uncertainty on the quantity.",
        verbose_name="Error low (-)",
    )

    error_sigma = models.FloatField(
        help_text="The confidence of the reported errors.",
        verbose_name="Error confidence",
        null=True,
        blank=True,
    )

    error_sigma_type = models.CharField(
        max_length=2,
        choices=[
            ('%', '%'),
            ('σ', 'σ'),
        ],
        default='σ',
        help_text="The type of the 'error_sigma' field",
        verbose_name="Error confidence type",
    )

    lower_limit = models.BooleanField(
        default=False,
        help_text="Is this an lower limit?",
    )

    upper_limit = models.BooleanField(
        default=False,
        help_text="Is this an upper limit?",
    )

    error_is_range = models.BooleanField(
        default=False,
        help_text="Display the error as a range, e.g. '1 - 5' instead of '3 ± 2'",
    )

    approximation = models.BooleanField(
        default=False,
        help_text="Display the error as a approximation, with the ~ symbol",
    )

    article = models.ForeignKey(
        "Article",
        help_text="The article in which this measurement is published",
        null=True,
        blank=True,
        related_name="measurements",
        on_delete=models.RESTRICT,
    )

    date = models.DateTimeField(
        null=True,
        blank=True,
        help_text="The date/time when this measurement was made",
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text="Any extra notes on this measurement",
    )

    owner = models.ForeignKey(
        User,
        null=True,
        blank=True,
        related_name="made_measurements",
        on_delete=models.RESTRICT,
        help_text="The user who made this measurement, if not a published measurement",
    )

    access = models.CharField(
        max_length=1,
        choices=[
            (ACCESS_PUBLIC, 'Public'),
            (ACCESS_GROUP, 'Group'),
            (ACCESS_PRIVATE, 'Only me'),
        ],
        default=ACCESS_PRIVATE,
        help_text="If published, everyone can access. If unpublished, this setting chooses who can see this measurement.",
    )

    access_groups = models.ManyToManyField(
        Group,
        blank=True,
        related_name="measurements",
        help_text="Which groups can view this measurement",
    )

    updated = models.DateTimeField(
        auto_now=True,
    )

    freq_ctr = models.FloatField(
        null=True,
        blank=True,
        help_text="The centre frequency.",
    )

    freq_hi = models.FloatField(
        null=True,
        blank=True,
        help_text="The top of the frequency range.",
    )

    freq_lo = models.FloatField(
        null=True,
        blank=True,
        help_text="The bottom of the frequency range.",
    )

    freq_band = models.ForeignKey(
        "FrequencyBand",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        help_text="The frequency band.",
        related_name="measurements",
    )

    freq_astropy_units = models.CharField(
        max_length=31,
        default="MHz",
        help_text="An astropy-conversant unit string that applies to the frequencies.",
    )

    chisq = models.FloatField(
        null=True,
        blank=True,
        verbose_name='χ²',
        help_text="The chi-squared of the fit.",
    )

    reduced_chisq = models.FloatField(
        null=True,
        blank=True,
        verbose_name='Reduced χ²',
        help_text="The reduced chi-squared of the fit.",
    )

    ANGLE_DDMMSS = 'D'
    ANGLE_HHMMSS = 'H'
    SPECTRAL_TYPE = 'S'

    special_display = models.CharField(
        max_length=1,
        null=True,
        blank=True,
        choices=[
            [ANGLE_DDMMSS, 'DDMMSS'],
            [ANGLE_HHMMSS, 'HHMMSS'],
            [SPECTRAL_TYPE, 'Harvard spectral classification'],
        ],
        help_text="Display in the chosen format.",
        verbose_name="Display format",
    )

    stokes = models.CharField(
        max_length=2,
        choices=[
            ['I', 'I'],
            ['Q', 'Q'],
            ['U', 'U'],
            ['V', 'V'],
        ],
        null=True,
        blank=True,
        help_text="The Stokes parameter associated with this measurement.",
    )

    notes = models.TextField(
        null=True,
        blank=True,
        help_text="Extra notes about this measurement.",
    )

    @property
    def astropy_quantity(self):
        return float(self.quantity) * 10**(self.power_of_10) * u.Unit(self.parameter.astropy_unit)

    @property
    def astropy_err(self):
        if self.err:
            return float(self.err) * 10**(self.power_of_10) * u.Unit(self.parameter.astropy_unit)

        return None

    @property
    def confidence(self):
        if not self.error_sigma:
            return None

        if self.error_sigma_type == '%':
            return f'{self.error_sigma:.0f}%'

        if self.error_sigma_type == 'σ':
            return f'{self.error_sigma:.0f}σ'

        return None


    @property
    def formatted_quantity(self):

        retstr = ""
        if self.precision is None:
            precision = Decimal(f"0.00000000000000000001")
        elif self.precision > 0:
            precision = Decimal(f"0.{'0'*((self.precision - 1) if self.precision else 19)}1")
        elif self.precision < 0:
            precision = Decimal(f"1{'0'*self.precision}")
        else:
            precision = Decimal("1")
        quantity = self.quantity.quantize(precision)

        if self.approximation:
            if self.upper_limit:
                retstr += "≲ "
            elif self.lower_limit:
                retstr += "≳ "
            else:
                retstr += "~ "
        else:
            if self.upper_limit:
                retstr += "< "
            elif self.lower_limit:
                retstr += "> "

        err    = self.err
        err_hi = self.err_hi
        err_lo = self.err_lo

        # Now that the prefix is taken care of, consider the special case
        # of spectral type. Must be dimensionless
        if self.parameter.astropy_unit is None:
            if self.special_display == self.SPECTRAL_TYPE:
                spectral_class = "OBAFGKM"[int(self.quantity/10)]
                spectral_subclass = self.quantity % 10
                retstr += spectral_class
                quantity = spectral_subclass.quantize(precision)

        if self.parameter.astropy_unit and u.Unit(self.parameter.astropy_unit).is_equivalent('deg'):
            if self.special_display == self.ANGLE_DDMMSS:
                angle = Angle(f'{self.quantity} {self.parameter.astropy_unit}').signed_dms
                retstr += f'{angle.d:02.0f}:{angle.m:02.0f}:'
                if round(angle.s) < 10: # zero padding for arcseconds
                    retstr += "0"
                quantity = Decimal(angle.s).quantize(precision)
                if err:
                    err = (self.err * 3600).quantize(precision)
                if err_hi:
                    err_hi = (self.err_hi * 3600).quantize(precision)
                if err_lo:
                    err_lo = (self.err_lo * 3600).quantize(precision)
            if self.special_display == self.ANGLE_HHMMSS:
                angle = Angle(f'{self.quantity} {self.parameter.astropy_unit}').hms
                retstr += f'{angle.h:02.0f}:{angle.m:02.0f}:'
                if round(angle.s) < 10: # zero padding for seconds
                    retstr += "0"
                quantity = Decimal(angle.s).quantize(precision)
                if err:
                    err = (self.err * 240).quantize(precision)
                if err_hi:
                    err_hi = (self.err_hi * 240).quantize(precision)
                if err_lo:
                    err_lo = (self.err_lo * 240).quantize(precision)

        if self.error_is_range == True and err is not None:
            lower_limit = (quantity - err).quantize(precision)
            upper_limit = (quantity + err).quantize(precision)

            quantity_str = f"{lower_limit} - {upper_limit}"

            if self.power_of_10 != 0:
                quantity_str = f"({quantity_str})"

        elif err:
            quantity_str = f"{quantity} ± {err.quantize(precision)}"

        else:
            quantity_str = f"{quantity}"
            if err_hi or err_lo:
                if err_hi == err_lo:
                    quantity_str = f" ± {err_hi.quantize(precision)}" # Basically, same as just "err" above
                else:
                    # CSS class supsub is defined in published/static/published/published-style.css
                    quantity_str += f"<span class='supsub'><sup>{err_hi.quantize(precision) if err_hi else '?'}</sup><sub>-{err_lo.quantize(precision) if err_lo else '?'}</sub></span>"

        if self.power_of_10 != 0:
            if err or err_hi or err_lo:
                quantity_str = f"({quantity_str})"
            if quantity != Decimal('1'):
                quantity_str += f" × 10<sup>{self.power_of_10}</sup>"
            else:
                quantity_str = f"10<sup>{self.power_of_10}</sup>"

        retstr += f"{quantity_str}"

        return retstr

    @property
    def formatted_quantity_with_units(self):

        if self.parameter.astropy_unit and u.Unit(self.parameter.astropy_unit).is_equivalent('deg'):
            if self.special_display == self.ANGLE_DDMMSS:
                return f"{self.formatted_quantity}″"
            if self.special_display == self.ANGLE_HHMMSS:
                return f"{self.formatted_quantity} s"

        if self.parameter.unicode_unit is not None:
            return f"{self.formatted_quantity} {self.parameter.unicode_unit}"

        return f"{self.formatted_quantity}"

    @property
    def accessible_by(self):
        if self.access == self.ACCESS_PUBLIC:
            return 'Everyone'

        if self.access == self.ACCESS_GROUP:
            return ', '.join([f'{access_group}' for access_group in self.access_groups.all()])

        # Else, only owner has access
        return f'{self.owner}'

    def check_access(self, user):
        # Return true if user has access to this measurement
        if self.access == self.ACCESS_PUBLIC:
            return True

        if self.access == self.ACCESS_GROUP:
            if self.access_groups.filter(user=user).exists():
                return True

        # Assert: self.access must = ACCESS_PRIVATE
        if self.user == user:
            return True

        # Default case
        return False

    @property
    def freq_band_display(self):
        #<td>{% if measurement.freq_ctr != None %}{{ measurement.freq_ctr }} {{ measurement.freq_astropy_units }}{% endif %}</td>
        if self.freq_band:
            return f'{self.freq_band.symbol}'
        if self.freq_lo and self.freq_hi:
            return f'{self.freq_lo} - {self.freq_hi} {self.freq_astropy_units}'
        if self.freq_ctr:
            return f'{self.freq_ctr} {self.freq_astropy_units}'
        return None

    def __str__(self):
        return self.formatted_quantity_with_units

    def __mul__(self, other):
        if instance(other, (int, float)):
            parameter = Parameter(
                astropy_unit=self.parameter.astropy_unit if self.parameter else None,
            )
            measurement = Measurement(
                quantity=other*self.quantity,
                parameter=parameter,
                power_of_10=self.power_of_10,
                err=(other*self.err) if self.err else None,
                err_hi=(other*self.err_hi) if self.err_hi else None,
                err_lo=(other*self.err_lo) if self.err_lo else None,
                err_sigma_type=self.error_sigma_type,
                lower_limit=self.lower_limit,
                upper_limit=self.upper_limit,
                error_is_range=self.error_is_range,
                approximation=self.approximation,
                special_display=self.special_display,
            )
            return measurement
        else:
            raise NotImplementedError

    __rmul__ = __mul__

    class Meta:
        ordering = ['quantity', 'ulp', 'parameter', 'article']


class ParameterSet(models.Model):

    name = models.CharField(
        max_length=63,
        unique=True,
        help_text="A name for the parameter set (e.g. 'ephemeris', 'main_table', 'viewing_geometry')",
    )

    description = models.TextField(
        blank=True,
        null=True,
        help_text="A helpful description of the parameter set",
    )

    parameters = models.ManyToManyField(
        "Parameter",
        blank=True,
        related_name="sets",
    )

    def __str__(self):
        return f"{self.name}"

    class Meta:
        ordering = ['name']

class Covariance(models.Model):

    measurement1 = models.ForeignKey(
        "Measurement",
        on_delete=models.CASCADE,
        help_text="The first measurement in the pair.",
        related_name="covariances_as_primary",
    )

    measurement2 = models.ForeignKey(
        "Measurement",
        on_delete=models.CASCADE,
        help_text="The second measurement in the pair.",
        related_name="covariances_as_secondary",
    )

    covariance = models.FloatField(
        help_text="The covariance value.",
    )

    def __str__(self):
        return f"{self.measurement1.parameter.name}, {self.measurement2.parameter.name}"

    class Meta:
        ordering = ['measurement1', 'measurement2']


class UserSetting(models.Model):

    SITE_THEMES = [
        ('l', 'Light'),
        ('d', 'Dark'),
    ]

    user = models.OneToOneField(
        User,
        on_delete=models.DO_NOTHING,
        related_name='setting',
    )

    site_theme = models.CharField(
        max_length=1,
        choices=SITE_THEMES,
        default='l',
        help_text="Choice of light or dark theme for the website",
    )

    def __str__(self) -> str:
        return f"Settings for {self.user}"


class Progenitor(models.Model):

    name = models.CharField(max_length=127, unique=True)
    abbr = models.CharField(max_length=31, verbose_name="Abbreviation", null=True, blank=True)

    def __str__(self) -> str:
        return self.abbr or self.name

    class Meta:
        ordering = ['name']


class ProgenitorClaim(models.Model):

    ulp = models.ForeignKey('Ulp', models.CASCADE, related_name='progenitor_claims')
    progenitor = models.ForeignKey('Progenitor', models.CASCADE, related_name='claims')
    article = models.ForeignKey('Article', models.CASCADE, related_name='progenitor_claims')
    tentative = models.IntegerField(default=0) # 0 = not tenatative, 1 = tentative (?), 2 = really tentative (??), etc.

    def __str__(self) -> str:
        return f"{self.progenitor}{'?'*self.tentative}"

    class Meta:
        ordering = ['ulp', 'progenitor']
        constraints = [
            models.UniqueConstraint(fields=['ulp', 'progenitor', 'article'], name='unique_ulp_progenitor_article')
        ]


class FrequencyBand(models.Model):

    FREQ_UNITS = u.GHz

    name = models.CharField(max_length=63, unique=True)
    symbol = models.CharField(max_length=15)
    lo = models.FloatField(help_text="Low end of frequency range (GHz).")
    hi = models.FloatField(help_text="High end of frequency range (GHz).")
    unit = models.CharField(max_length=15, help_text="Display unit. AstroPy.unit-conformant unit in units of either length (for wavelength) or inverse time (for frequency).")

    def __str__(self) -> str:
        return f'{self.name}'

    @property
    def range_display(self):
        f1 = (self.lo*self.FREQ_UNITS).to(self.unit, equivalencies=u.spectral())
        f2 = (self.hi*self.FREQ_UNITS).to(self.unit, equivalencies=u.spectral())
        lo, hi = (f1, f2) if f1 < f2 else (f2, f1)
        return f'{lo.value}-{hi}'

    class Meta:
        ordering = ['lo', 'hi']
