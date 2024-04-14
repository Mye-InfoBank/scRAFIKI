FROM python:3.9.18-slim-bullseye

# Install git
RUN apt-get update && apt-get install -y git build-essential procps graphviz

RUN pip install --no-cache-dir --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html && \
    pip install --no-cache-dir scarches \
        git+https://github.com/Mye-InfoBank/scib.git#egg=scib[louvain,desc,bbknn,scanorama,scgen,scvi,trvaep,harmony]
