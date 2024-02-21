# Generated by Django 5.0.2 on 2024-02-21 12:39

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):
    initial = True

    dependencies = [
        ("auth", "0012_alter_user_first_name_max_length"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name="Article",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                (
                    "citet_text",
                    models.CharField(
                        blank=True,
                        help_text="The display text when this is cited inline, e.g., Smith et. al (1900)",
                        max_length=127,
                        null=True,
                    ),
                ),
                (
                    "title",
                    models.CharField(
                        blank=True,
                        help_text="The title of the paper",
                        max_length=1023,
                        null=True,
                    ),
                ),
                (
                    "doi",
                    models.CharField(
                        blank=True,
                        help_text="The DOI of the article",
                        max_length=1023,
                        null=True,
                        verbose_name="DOI",
                    ),
                ),
            ],
            options={
                "ordering": ["citet_text"],
            },
        ),
        migrations.CreateModel(
            name="Parameter",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                (
                    "name",
                    models.CharField(
                        help_text="The name of the parameter",
                        max_length=63,
                        unique=True,
                    ),
                ),
                (
                    "description",
                    models.TextField(
                        blank=True,
                        help_text="A description of this measurable parameter",
                        null=True,
                    ),
                ),
                (
                    "ascii_symbol",
                    models.CharField(
                        help_text="The symbol used for this parameter, expressed in ASCII characters",
                        max_length=15,
                    ),
                ),
                (
                    "latex_symbol",
                    models.CharField(
                        blank=True,
                        help_text="The symbol used for this parameter, expressed as a LaTeX expression",
                        max_length=63,
                        null=True,
                    ),
                ),
                (
                    "astropy_unit",
                    models.CharField(
                        blank=True,
                        help_text="The (astropy-conversant) unit of this parameter (e.g. 's', 'lightyear')",
                        max_length=31,
                        null=True,
                    ),
                ),
            ],
            options={
                "ordering": ["name"],
            },
        ),
        migrations.CreateModel(
            name="Measurement",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                (
                    "quantity",
                    models.DecimalField(
                        decimal_places=20,
                        help_text="The numerical value of this measurement",
                        max_digits=40,
                    ),
                ),
                ("power_of_10", models.IntegerField(default=0, verbose_name="× 10^")),
                (
                    "precision",
                    models.IntegerField(
                        blank=True,
                        help_text="The number of decimal places to display",
                        null=True,
                    ),
                ),
                (
                    "err",
                    models.DecimalField(
                        blank=True,
                        decimal_places=20,
                        help_text="The (symmetrical) uncertainty on the quantity.",
                        max_digits=40,
                        null=True,
                        verbose_name="Error (±)",
                    ),
                ),
                (
                    "err_hi",
                    models.DecimalField(
                        blank=True,
                        decimal_places=20,
                        help_text="The (asymmetrical) upper uncertainty on the quantity.",
                        max_digits=40,
                        null=True,
                        verbose_name="Error high (+)",
                    ),
                ),
                (
                    "err_lo",
                    models.DecimalField(
                        blank=True,
                        decimal_places=20,
                        help_text="The (asymmetrical) lower uncertainty on the quantity.",
                        max_digits=40,
                        null=True,
                        verbose_name="Error low (-)",
                    ),
                ),
                (
                    "error_sigma",
                    models.FloatField(
                        default=1,
                        help_text="The sigma of the reported errors.",
                        verbose_name="Error significance, σ",
                    ),
                ),
                (
                    "lower_limit",
                    models.BooleanField(
                        default=False, help_text="Is this an lower limit?"
                    ),
                ),
                (
                    "upper_limit",
                    models.BooleanField(
                        default=False, help_text="Is this an upper limit?"
                    ),
                ),
                (
                    "error_is_range",
                    models.BooleanField(
                        default=False,
                        help_text="Display the error as a range, e.g. '1 - 5' instead of '3 ± 2'",
                    ),
                ),
                (
                    "approximation",
                    models.BooleanField(
                        default=False,
                        help_text="Display the error as a approximation, with the ~ symbol",
                    ),
                ),
                (
                    "date",
                    models.DateTimeField(
                        blank=True,
                        help_text="The date/time when this measurement was made",
                        null=True,
                    ),
                ),
                (
                    "access",
                    models.CharField(
                        choices=[("P", "Public"), ("G", "Group"), ("O", "Only me")],
                        default="O",
                        help_text="If published, everyone can access. If unpublished, this setting chooses who can see this measurement.",
                        max_length=1,
                    ),
                ),
                ("updated", models.DateTimeField(auto_now=True)),
                (
                    "freq_ctr",
                    models.FloatField(
                        blank=True, help_text="The centre frequency.", null=True
                    ),
                ),
                (
                    "freq_hi",
                    models.FloatField(
                        blank=True,
                        help_text="The top of the frequency range.",
                        null=True,
                    ),
                ),
                (
                    "freq_lo",
                    models.FloatField(
                        blank=True,
                        help_text="The bottom of the frequency range.",
                        null=True,
                    ),
                ),
                (
                    "freq_astropy_units",
                    models.CharField(
                        default="MHz",
                        help_text="An astropy-conversant unit string that applies to the frequencies.",
                        max_length=31,
                    ),
                ),
                (
                    "chisq",
                    models.FloatField(
                        blank=True,
                        help_text="The chi-squared of the fit.",
                        null=True,
                        verbose_name="χ²",
                    ),
                ),
                (
                    "reduced_chisq",
                    models.FloatField(
                        blank=True,
                        help_text="The reduced chi-squared of the fit.",
                        null=True,
                        verbose_name="Reduced χ²",
                    ),
                ),
                (
                    "angle_display",
                    models.CharField(
                        blank=True,
                        choices=[("D", "DDMMSS"), ("H", "HHMMSS")],
                        help_text="If an angle, display in the chosen sexagesimal format.",
                        max_length=1,
                        null=True,
                    ),
                ),
                (
                    "stokes",
                    models.CharField(
                        blank=True,
                        choices=[("I", "I"), ("Q", "Q"), ("U", "U"), ("V", "V")],
                        help_text="The Stokes parameter associated with this measurement.",
                        max_length=2,
                        null=True,
                    ),
                ),
                (
                    "notes",
                    models.TextField(
                        blank=True,
                        help_text="Extra notes about this measurement.",
                        null=True,
                    ),
                ),
                (
                    "access_groups",
                    models.ManyToManyField(
                        blank=True,
                        help_text="Which groups can view this measurement",
                        related_name="measurements",
                        to="auth.group",
                    ),
                ),
                (
                    "article",
                    models.ForeignKey(
                        blank=True,
                        help_text="The article in which this measurement is published",
                        null=True,
                        on_delete=django.db.models.deletion.RESTRICT,
                        related_name="measurements",
                        to="published.article",
                    ),
                ),
                (
                    "owner",
                    models.ForeignKey(
                        blank=True,
                        help_text="The user who made this measurement, if not a published measurement",
                        null=True,
                        on_delete=django.db.models.deletion.RESTRICT,
                        related_name="made_measurements",
                        to=settings.AUTH_USER_MODEL,
                    ),
                ),
                (
                    "parameter",
                    models.ForeignKey(
                        help_text="The parameter beaing measured",
                        on_delete=django.db.models.deletion.RESTRICT,
                        related_name="measurements",
                        to="published.parameter",
                    ),
                ),
            ],
            options={
                "ordering": ["ulp", "article", "parameter"],
            },
        ),
        migrations.CreateModel(
            name="Covariance",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("covariance", models.FloatField(help_text="The covariance value.")),
                (
                    "measurement1",
                    models.ForeignKey(
                        help_text="The first measurement in the pair.",
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="covariances_as_primary",
                        to="published.measurement",
                    ),
                ),
                (
                    "measurement2",
                    models.ForeignKey(
                        help_text="The second measurement in the pair.",
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="covariances_as_secondary",
                        to="published.measurement",
                    ),
                ),
            ],
            options={
                "ordering": ["measurement1", "measurement2"],
            },
        ),
        migrations.CreateModel(
            name="ParameterSet",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                (
                    "name",
                    models.CharField(
                        help_text="A name for the parameter set (e.g. 'ephemeris', 'main_table', 'viewing_geometry')",
                        max_length=63,
                        unique=True,
                    ),
                ),
                (
                    "description",
                    models.TextField(
                        blank=True,
                        help_text="A helpful description of the parameter set",
                        null=True,
                    ),
                ),
                (
                    "parameters",
                    models.ManyToManyField(
                        blank=True, related_name="sets", to="published.parameter"
                    ),
                ),
            ],
            options={
                "ordering": ["name"],
            },
        ),
        migrations.CreateModel(
            name="Ulp",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                (
                    "name",
                    models.CharField(
                        help_text="The official name of the ULP",
                        max_length=63,
                        unique=True,
                    ),
                ),
                (
                    "abbr",
                    models.CharField(
                        blank=True,
                        help_text="A short version of the name",
                        max_length=15,
                        null=True,
                        unique=True,
                    ),
                ),
                (
                    "discoverer",
                    models.ForeignKey(
                        blank=True,
                        help_text="The group who first discovered this ULP",
                        null=True,
                        on_delete=django.db.models.deletion.RESTRICT,
                        related_name="ulps",
                        to="auth.group",
                    ),
                ),
            ],
            options={
                "verbose_name": "ULP",
                "verbose_name_plural": "ULPs",
                "ordering": ["name"],
            },
        ),
        migrations.AddField(
            model_name="measurement",
            name="ulp",
            field=models.ForeignKey(
                help_text="The ULP for to which this measurement applies",
                on_delete=django.db.models.deletion.RESTRICT,
                related_name="measurements",
                to="published.ulp",
                verbose_name="ULP",
            ),
        ),
    ]
