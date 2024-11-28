# Generated by Django 5.0.8 on 2024-11-28 00:59

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0043_toa_baseline_level_toa_baseline_slope'),
    ]

    operations = [
        migrations.AddField(
            model_name='toa',
            name='include_in_fit',
            field=models.BooleanField(default=True, help_text='Whether to include this ToA in the timing fit.'),
        ),
    ]
