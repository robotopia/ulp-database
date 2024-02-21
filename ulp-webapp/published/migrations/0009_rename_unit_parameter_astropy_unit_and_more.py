# Generated by Django 4.2.6 on 2024-02-09 02:16

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("published", "0008_parameterset_description"),
    ]

    operations = [
        migrations.RenameField(
            model_name="parameter",
            old_name="unit",
            new_name="astropy_unit",
        ),
        migrations.AlterField(
            model_name="parameterset",
            name="name",
            field=models.CharField(
                help_text="A name for the parameter set (e.g. 'ephemeris', 'main_table', 'viewing_geometry')",
                max_length=63,
                unique=True,
            ),
        ),
    ]