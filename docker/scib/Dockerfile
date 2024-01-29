FROM python:3.9.18-slim-bullseye

ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR
ENV CELLTYPIST_FOLDER=/tmp/celltypist
RUN mkdir -p $CELLTYPIST_FOLDER && chmod -R 777 $CELLTYPIST_FOLDER
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p $MPLCONFIGDIR && chmod -R 777 $MPLCONFIGDIR

# Install git
RUN apt-get update && apt-get install -y git build-essential procps graphviz

RUN pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

RUN pip install --no-cache-dir \
        git+https://github.com/Mye-InfoBank/scib.git#egg=scib[louvain,desc,bbknn,scanorama,scgen,scvi,trvaep,harmony]
