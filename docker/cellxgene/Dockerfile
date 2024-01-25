FROM ubuntu:latest

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update && \
    apt-get install -y build-essential libxml2-dev python3-dev python3-pip zlib1g-dev python3-requests python3-aiohttp && \
    python3 -m pip install --upgrade pip && \
    pip3 install cellxgene 

WORKDIR /srv

ENTRYPOINT ["cellxgene"]

# use this command to start the container without docker-compose.yml. The paths and ports need to be adapted to your local system.
#docker run --rm --network='host' -p 5456:5456 -v ~/cost/atlas:/atlas cellxgene launch /atlas/adata_cell_type_coarse.h5ad -p 5456