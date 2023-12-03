FROM rocker/r2u:22.04

ENV WORKON_HOME=/opt/venv

RUN apt-get update && apt install -y git libffi-dev libssl-dev libsqlite3-dev curl python3

COPY install.R /tmp/
RUN Rscript /tmp/install.R