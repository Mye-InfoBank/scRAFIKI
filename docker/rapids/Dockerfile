FROM mambaorg/micromamba:1.5-jammy

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/bin:$PATH