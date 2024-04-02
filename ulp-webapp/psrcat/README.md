# Importing ATNF pulsar catalogue v2.1 into a Django webapp

1. Download the database (`psrcat2.db`) from [CSIRO](https://doi.org/10.25919/d9a8-pq84) into the webapp's project folder (`ulp-webapp/ulp-webapp`).
2. Add an extra entry to `DATABASES` in `ulp-webapp/settings.py`:
```
DATABASES = {
    ...
    'psrcat': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(PROJECT_DIR, 'psrcat2.db'),
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    }
}
```
Also, add `psrcat` to `INSTALLED_APPS`:
```
    'psrcat.apps.PsrcatConfig',
```

3. From the parent directory (`ulp-webapp`), use `inspectdb` to create models:
```
python manage.py inspectdb --database=psrcat > psrcat/models.py
```
From the generated `models.py` file, manually go through and remove the `blank=True, null=True` options for all fields marked with `primary_key=True`.
Optionally, change the `TextField` fields to `CharField`, remembering to add an appropriate `max_length` argument.
