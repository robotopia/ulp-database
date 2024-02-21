#!/bin/bash
set -e

# Initialise the DB with the following 
# Have to do this per line otherwise postgres complains
psql -v ON_ERROR_STOP=1 -c "CREATE DATABASE $DBNAME;"
psql -v ON_ERROR_STOP=1 -c "CREATE USER $DBUSER WITH ENCRYPTED PASSWORD '$DBPASS';"
psql -v ON_ERROR_STOP=1 -c "ALTER ROLE $DBUSER SET client_encoding TO 'utf8';"
psql -v ON_ERROR_STOP=1 -c "ALTER ROLE $DBUSER SET default_transaction_isolation TO 'read committed';"
psql -v ON_ERROR_STOP=1 -c "ALTER ROLE $DBUSER SET timezone TO 'UTC';"
