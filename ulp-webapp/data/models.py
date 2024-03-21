from django.db import models
from django.contrib.auth.models import User, Group
import astropy.units as u
from astropy.coordinates import Angle
from decimal import Decimal

# Create your models here.
class Telescope(models.Model):

    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="The name of the telescope",
    )

    abbr = models.CharField(
        max_length=15,
        unique=True,
        blank=True,
        null=True,
        help_text="The abbreviation of the name",
        verbose_name="Abbreviation",
    )

    latitude_deg = models.DecimalField(
        decimal_places=10,
        max_digits=14,
        null=True,
        blank=True,
        help_text="The latitude (in degrees) of the telescope",
        verbose_name="Latitude (°)",
    )

    longitude_deg = models.DecimalField(
        decimal_places=10,
        max_digits=14,
        null=True,
        blank=True,
        help_text="The longitude (in degrees) of the telescope",
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
        help_text="The telescope to which this backend belongs",
        related_name="backends",
    )

    name = models.CharField(
        max_length=255,
        unique=True,
        help_text="The name of the telescope backend",
    )

    freq_ctr = models.FloatField(
        help_text="The centre frequency of the band",
        verbose_name="Centre frequency",
    )

    bw = models.FloatField(
        help_text="The bandwidth",
        verbose_name="Bandwidth",
    )

    freq_units = models.CharField(
        max_length=31,
        default="MHz",
        help_text="An astropy-conversant unit string that applies to the frequencies.",
    )

    def __str__(self):
        return f'{self.name} ({self.telescope})'

    class Meta:
        ordering = ['telescope', 'name',]

