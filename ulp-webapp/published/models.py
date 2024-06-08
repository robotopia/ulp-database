from django.db import models
from django.contrib.auth.models import User, Group

from decimal import Decimal

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord

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
        help_text="The display text when this is cited inline, e.g., Smith et. al (1900)",
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

        return u.Unit(self.astropy_unit).to_string(format='unicode')

    class Meta:
        ordering = ['name']


class Measurement(models.Model):

    ACCESS_PUBLIC = 'P'
    ACCESS_GROUP = 'G'
    ACCESS_PRIVATE = 'O'

    parameter = models.ForeignKey(
        "Parameter",
        help_text="The parameter beaing measured",
        related_name="measurements",
        on_delete=models.RESTRICT,
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
        default=1,
        help_text="The sigma of the reported errors.",
        verbose_name="Error significance, σ",
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

    angle_display = models.CharField(
        max_length=1,
        null=True,
        blank=True,
        choices=[
            [ANGLE_DDMMSS, 'DDMMSS'],
            [ANGLE_HHMMSS, 'HHMMSS'],
        ],
        help_text="If an angle, display in the chosen sexagesimal format.",
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
    def formatted_quantity(self):

        if self.parameter.astropy_unit and u.Unit(self.parameter.astropy_unit).is_equivalent('deg'):
            if self.angle_display == self.ANGLE_DDMMSS:
                return Angle(f'{self.quantity} {self.parameter.astropy_unit}').to_string(unit=u.deg, pad=True, format='unicode')
            if self.angle_display == self.ANGLE_HHMMSS:
                return Angle(f'{self.quantity} {self.parameter.astropy_unit}').to_string(unit=u.hourangle, pad=True, format='unicode')

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

        if self.error_is_range == True and self.err is not None:
            lower_limit = (self.quantity - self.err).quantize(precision)
            upper_limit = (self.quantity + self.err).quantize(precision)

            quantity_str = f"{lower_limit} - {upper_limit}"

            if self.power_of_10 != 0:
                quantity_str = f"({quantity_str})"

        elif self.err:
            err = self.err.quantize(precision)
            quantity_str = f"{quantity} ± {err}"

        else:
            quantity_str = f"{quantity}"
            if self.err_hi:
                err_hi = self.err_hi.quantize(precision)
                quantity_str += f"^+{err_hi}"
            if self.err_lo:
                err_lo = self.err_lo.quantize(precision)
                quantity_str += f"^+{err_lo}"

        if self.power_of_10 != 0:
            if self.err or self.err_hi or self.err_lo:
                quantity_str = f"({quantity_str})"
            if quantity != Decimal('1'):
                quantity_str += f" × 10^{self.power_of_10}"
            else:
                quantity_str = f"10^{self.power_of_10}"

        retstr += f"{quantity_str}"

        return retstr

    @property
    def formatted_quantity_with_units(self):

        if self.parameter.astropy_unit and u.Unit(self.parameter.astropy_unit).is_equivalent('deg'):
            return f"{self.formatted_quantity}"

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

    def __str__(self):
        return self.formatted_quantity_with_units

    class Meta:
        ordering = ['ulp', 'article', 'parameter']


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
        return f"{measurement1.parameter.name}, {measurement2.parameter.name}"

    class Meta:
        ordering = ['measurement1', 'measurement2']
