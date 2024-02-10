# Generated by Django 5.0.2 on 2024-02-10 00:35

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        (
            "published",
            "0013_measurement_freq_astropy_units_measurement_freq_ctr_and_more",
        ),
    ]

    operations = [
        migrations.AlterModelOptions(
            name="article",
            options={"ordering": ["citet_text"]},
        ),
        migrations.AlterModelOptions(
            name="parameter",
            options={"ordering": ["name"]},
        ),
        migrations.AlterModelOptions(
            name="ulp",
            options={
                "ordering": ["name"],
                "verbose_name": "ULP",
                "verbose_name_plural": "ULPs",
            },
        ),
        migrations.RemoveField(
            model_name="measurement",
            name="prefer_scientific_notation",
        ),
        migrations.RemoveField(
            model_name="parameter",
            name="prefer_scientific_notation",
        ),
    ]
