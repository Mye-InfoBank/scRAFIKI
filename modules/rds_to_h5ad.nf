process RDS_TO_H5AD {
  tag "${meta.id}"

  container "bigdatainbiomedicine/sc-rpy:1.0"
  label "process_medium"

  input:
  tuple val(meta), path(rds)
  
  output:
  tuple val(meta), path("${meta.id}.h5ad")
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import anndata2ri
  import rpy2.robjects as ro
  import numpy as np
  seurat = ro.packages.importr('Seurat')

  sce = ro.r(f'as.SingleCellExperiment(readRDS("${rds}"))')

  adata = anndata2ri.rpy2py(sce)

  # Convert indices to string
  adata.obs.index = adata.obs.index.astype(str)
  adata.var.index = adata.var.index.astype(str)

  adata.write_h5ad("${meta.id}.h5ad")
  """
}