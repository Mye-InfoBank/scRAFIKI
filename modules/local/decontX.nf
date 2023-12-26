process DECONTX {
  tag "$meta.id"

  container "bigdatainbiomedicine/sc-r"

  label "process_medium"

  input:
    tuple val(meta), file(adata)
  
  output:
    tuple val(meta), file("*_ambient.h5ad")
  
  script:
  """
  #!/usr/bin/env Rscript

  library(anndata)
  library(celda)

  file_path <- "$adata"

  ad <- read_h5ad(file_path)

  results <- decontX(t(ad\$X))

  ad\$X <- t(results\$decontXcounts)

  basename <- gsub(".h5ad", "", basename(file_path))

  write_h5ad(ad, paste0(basename, "_ambient.h5ad"))
  """
}