FROM mambaorg/micromamba:1.5-jammy

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.9 git graphviz celltypist && \
    micromamba run -n base pip install --no-cache-dir \
        git+https://github.com/Mye-InfoBank/scib.git#egg=scib[louvain,desc,bbknn,scanorama,scgen,scvi,trvaep,harmony] && \
    micromamba clean --all --yes

ENV PATH /opt/conda/bin:$PATH

WORKDIR /app

ENV NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p $NUMBA_CACHE_DIR && chmod -R 777 $NUMBA_CACHE_DIR
ENV CELLTYPIST_FOLDER=/tmp/celltypist
RUN mkdir -p $CELLTYPIST_FOLDER && chmod -R 777 $CELLTYPIST_FOLDER
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p $MPLCONFIGDIR && chmod -R 777 $MPLCONFIGDIR