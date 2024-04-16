FROM continuumio/miniconda3:24.1.2-0

COPY environment.yml /tmp/environment.yml

RUN conda env update -n base -f /tmp/environment.yml && \
    conda clean -afy