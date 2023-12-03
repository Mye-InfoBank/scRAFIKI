install.packages(c("BiocManager", "Seurat", "anndata", "Matrix"))
BiocManager::install(c("SingleR", "decontX"))

reticulate::install_python()

reticulate::virtualenv_create()

anndata::install_anndata()