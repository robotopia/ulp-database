from django.db import models

# Create your models here.
class ulp(models.Model):

    name = models.CharField(
        max_length=63,
        unique=True,
        help_text="The official name of the ULP",
    )

    class Meta:
        verbose_name = "ULP"
        verbose_name_plural = "ULPs"
