FROM python:3.11

# Copy webapp into container
RUN mkdir /ulp-webapp
WORKDIR /ulp-webapp
COPY . /ulp-webapp/

# Install necessary packages
RUN apt-get update
RUN apt-get install build-essential f2c postgresql libpq-dev python3.11-dev python3-pip gosu locales vim -yqq

# Install Python requirements
RUN pip3 install uwsgi
RUN pip3 install -r requirements.txt

# Start the container with entrypoint script
ENTRYPOINT ["/ulp-webapp/entrypoint.sh"]
