FROM rocker/r2u:latest

# Prepare numba cache directory
ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR

WORKDIR /tmp

# Install python 3.9
RUN apt-get update && apt-get install -y wget openssl libssl-dev libffi-dev libsqlite3-dev && \
    wget https://www.python.org/ftp/python/3.9.0/Python-3.9.0.tgz && \
    tar -xvf Python-3.9.0.tgz && \
    cd Python-3.9.0 && \
    ./configure --enable-optimizations && \
    make -j 4 && \
    make install

# Install python and R packages for pyScTransform
COPY requirements.txt install.R /tmp/

RUN Rscript install.R
RUN python3 -m pip install -r requirements.txt

WORKDIR /app