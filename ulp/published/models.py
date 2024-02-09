from django.db import models
from django.contrib.auth.models import User, Group
import astropy.units as u
from decimal import Decimal

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

    readers = models.ManyToManyField(
        User,
        blank=True,
        help_text="The users that have read access to this ULP",
        related_name="readable_ulps",
    )

    editors = models.ManyToManyField(
        User,
        blank=True,
        help_text="The users that have edit access to this ULP",
        related_name="editable_ulps",
    )

    class Meta:
        verbose_name = "ULP"
        verbose_name_plural = "ULPs"

    def __str__(self):
        return self.name


class Article(models.Model):

    citet_text = models.CharField(
        max_length=127,
        blank=True,
        null=True,
        help_text="The display text when this is cited inline, e.g., Smith et. al (1900)",
    )

    doi = models.CharField(
        max_length=1023,
        blank=True,
        null=True,
        help_text="The DOI of the article",
        verbose_name="DOI",
    )

    def __str__(self):
        return self.citet_text


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

    prefer_scientific_notation = models.BooleanField(
        default=False,
        help_text="True: 3x10²; False: 300",
    )

    def __str__(self):
        return f"{self.name} ({u.Unit(self.astropy_unit).to_string(format='unicode')})"


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
        help_text="The true quantity is multiplied by this power of 10. Same with the errors.",
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
        help_text="The (symmetrical) uncertainty on the quantity",
    )

    err_hi = models.DecimalField(
        decimal_places=20,
        max_digits=40,
        null=True,
        blank=True,
        help_text="The (asymmetrical) upper uncertainty on the quantity",
    )

    err_lo = models.DecimalField(
        decimal_places=20,
        max_digits=40,
        null=True,
        blank=True,
        help_text="The (asymmetrical) lower uncertainty on the quantity",
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

    prefer_scientific_notation = models.BooleanField(
        default=False,
        help_text="True: 3x10²; False: 300",
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
            quantity_str = f"{quantity_str} × 10^{self.power_of_10}"

        retstr += f"{quantity_str}"

        return retstr

    @property
    def formatted_quantity_with_units(self):

        return f"{self.formatted_quantity} {self.parameter.astropy_unit}"

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
