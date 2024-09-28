# Generated by Django 5.0.8 on 2024-09-28 08:10

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("data", "0016_remove_timeofarrival_freq_units_and_more"),
    ]

    operations = [
        migrations.AlterField(
            model_name="timeofarrival",
            name="freq",
            field=models.FloatField(
                blank=True,
                help_text="The centre frequency of this detection, in MHz.",
                null=True,
                verbose_name="Frequency (MHz)",
            ),
        ),
    ]
