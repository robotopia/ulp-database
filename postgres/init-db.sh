#!/bin/bash
set -e

# Initialise the DB with the following 
psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB"  <<-EOSQL
    CREATE DATABASE $DBNAME;
    CREATE USER $DBUSER WITH ENCRYPTED PASSWORD '$DBPASS';
    ALTER ROLE $DBUSER SET client_encoding TO 'utf8';
    ALTER ROLE $DBUSER SET default_transaction_isolation TO 'read committed';
    ALTER ROLE $DBUSER SET timezone TO 'UTC';
    GRANT ALL PRIVILEGES ON DATABASE $DBNAME TO $DBUSER;
EOSQL
