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

    psql "postgres://${DBUSER}:${DBPASS}@postgres:5432/${DBNAME}" -v ON_ERROR_STOP=1 -c "\
        CREATE EXTENSION q3c; \
        CREATE INDEX ON candidate_app_candidate (q3c_ang2ipix(ra_deg, dec_deg)); \
    "
fi 

# This runs the web app locally through Django
# python3 manage.py runserver 0.0.0.0:80 

# This runs the webapp using uwsgi and creates a socket that nginx uses
uwsgi --ini /ulp-webapp/ulp-webapp.uwsgi.ini
