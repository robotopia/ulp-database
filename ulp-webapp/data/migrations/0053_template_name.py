# Generated by Django 5.0.10 on 2024-12-19 02:13

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0052_workingephemeris_tausc_1ghz'),
    ]

    operations = [
        migrations.AddField(
            model_name='template',
            name='name',
            field=models.CharField(blank=True, help_text='An optional name to help the user distinguish multiple templates for the same source.', max_length=127, null=True),
        ),
    ]
