# Generated by Django 5.0.8 on 2024-10-02 06:48

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0027_pulse'),
    ]

    operations = [
        migrations.AddField(
            model_name='timeofarrival',
            name='pulse',
            field=models.ForeignKey(blank=True, help_text='The pulse from which this ToA was derived.', null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='toas', to='data.pulse'),
        ),
    ]
