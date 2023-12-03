FROM python:3.11

COPY requirements.txt setup_env.sh /tmp/
COPY lib /tmp/lib

WORKDIR /tmp

# Install dependencies
RUN bash setup_env.sh

WORKDIR /app

# Remove the temporary files
RUN rm -rf /tmp

ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR
ENV CELLTYPIST_FOLDER=/tmp/celltypist
RUN mkdir -p $CELLTYPIST_FOLDER && chmod -R 777 $CELLTYPIST_FOLDER