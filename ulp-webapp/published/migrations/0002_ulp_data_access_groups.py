# Generated by Django 5.0.2 on 2024-03-21 08:02

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('auth', '0012_alter_user_first_name_max_length'),
        ('published', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='ulp',
            name='data_access_groups',
            field=models.ManyToManyField(blank=True, help_text="The groups that are allowed access to this ULP's data.", related_name='data_accessible_ulps', to='auth.group'),
        ),
    ]
