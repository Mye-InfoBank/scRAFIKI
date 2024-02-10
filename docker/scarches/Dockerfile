FROM python:3.9.18-slim-bullseye

ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p $MPLCONFIGDIR && chmod -R 777 $MPLCONFIGDIR

RUN apt-get update && apt-get install -y procps

RUN pip install scarches

RUN pip install --upgrade "jax[cuda12_pip]==0.4.20" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
