#!/bin/bash
set -e

# This is only run on the first start-up of the container.
CONTAINER_FIRST_STARTUP="CONTAINER_FIRST_STARTUP"
if [ ! -e /$CONTAINER_FIRST_STARTUP ]; then
    touch /$CONTAINER_FIRST_STARTUP

    # Iniitalise the DB with schema
    python3 manage.py makemigrations published
    python3 manage.py migrate published
    python3 manage.py migrate
    python3 manage.py migrate --run-syncdb
    python3 manage.py collectstatic

    # Download a planetary ephemeris to this folder
    wget -O de430.bsp http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp

fi 

if [ "$DJANGO_DEBUG" == "True" ]
then
    # This runs the web app locally through Django
    python3 manage.py runserver 0.0.0.0:8000
else
    # This runs the webapp using uwsgi and creates a socket that nginx uses
    uwsgi --ini /ulp-webapp/ulp-webapp.uwsgi.ini
fi
