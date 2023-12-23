FROM mambaorg/micromamba:1.5-jammy

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python \
        pandas anndata scanpy \
        louvain leidenalg keras \
        mnnpy harmony-pytorch scvi-tools scib \
        seaborn scikit-misc celltypist scanorama bbknn && \
    micromamba run -n base pip install --no-cache-dir \
        desc trvae trvaep scgen && \
    micromamba clean --all --yes

ENV PATH /opt/conda/bin:$PATH

WORKDIR /app

ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR
ENV CELLTYPIST_FOLDER=/tmp/celltypist
RUN mkdir -p $CELLTYPIST_FOLDER && chmod -R 777 $CELLTYPIST_FOLDER
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p $MPLCONFIGDIR && chmod -R 777 $MPLCONFIGDIR