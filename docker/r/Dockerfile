FROM rocker/r2u:latest

RUN Rscript -e 'install.packages(c("BiocManager", "Seurat", "anndata", "MatrixExtra"))'

RUN Rscript -e 'BiocManager::install(c("celda", "SingleR"))'

