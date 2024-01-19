#!/bin/bash

conda install -c anaconda -c rapidsai -c bioconda -c conda-forge -c nvidia \
    rapids=23.12 python=3.10 cuda-version=12.0 wget scikit-misc \
    omnipath gdown cudnn cutensor cusparselt leidenalg louvain \
    fa2 pip harmonypy