# Generated by Django 5.0.2 on 2024-02-08 09:49

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('auth', '0012_alter_user_first_name_max_length'),
        ('published', '0002_article_parameter_ulp_abbr_measurement'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.AddField(
            model_name='ulp',
            name='discoverer',
            field=models.ForeignKey(blank=True, help_text='The group who first discovered this ULP', null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='ulps', to='auth.group'),
        ),
        migrations.AddField(
            model_name='ulp',
            name='editors',
            field=models.ManyToManyField(blank=True, help_text='The users that have edit access to this ULP', related_name='editable_ulps', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='ulp',
            name='readers',
            field=models.ManyToManyField(blank=True, help_text='The users that have read access to this ULP', related_name='readable_ulps', to=settings.AUTH_USER_MODEL),
        ),
    ]