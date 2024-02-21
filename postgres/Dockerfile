FROM postgres:14

# This init script is only run once at first run of the container.
COPY ./init-db.sh /docker-entrypoint-initdb.d/init-db.sh
