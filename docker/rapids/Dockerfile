FROM mambaorg/micromamba:1.5-jammy

ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR
ENV CELLTYPIST_FOLDER=/tmp/celltypist
RUN mkdir -p $CELLTYPIST_FOLDER && chmod -R 777 $CELLTYPIST_FOLDER
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p $MPLCONFIGDIR && chmod -R 777 $MPLCONFIGDIR
ENV CUPY_CACHE_DIR=/tmp/cupy_cache
RUN mkdir -p $CUPY_CACHE_DIR && chmod -R 777 $CUPY_CACHE_DIR

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /tmp/requirements.txt
RUN micromamba run -n base pip install -r /tmp/requirements.txt && \
    micromamba clean --all --yes

ENV PATH /opt/conda/bin:$PATH
