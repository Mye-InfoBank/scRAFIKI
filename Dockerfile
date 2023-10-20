FROM python:3.11

COPY requirements.txt setup_env.sh /tmp/
COPY lib /tmp/lib

WORKDIR /tmp

# Install dependencies
RUN bash setup_env.sh

WORKDIR /app
