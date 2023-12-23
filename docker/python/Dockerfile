FROM ubuntu:20.04

# manually enable gcc8 installation 
RUN echo 'deb [arch=amd64] http://archive.ubuntu.com/ubuntu focal main universe' >> /etc/apt/sources.list

# make sure gcc, python3 (dev) and pip are up to date
RUN apt-get update && apt-get install -y \
    gcc-8 g++-8 python3 python3-pip \ 
    python3-dev 

COPY requirements.txt /tmp/

WORKDIR /tmp

# Install dependencies
RUN pip install -r requirements.txt

WORKDIR /app

# Remove the temporary files
RUN rm -rf /tmp

ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR
ENV CELLTYPIST_FOLDER=/tmp/celltypist
RUN mkdir -p $CELLTYPIST_FOLDER && chmod -R 777 $CELLTYPIST_FOLDER
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p $MPLCONFIGDIR && chmod -R 777 $MPLCONFIGDIR