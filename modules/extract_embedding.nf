process EXTRACT_EMBEDDING {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  tuple val(meta), path(input)
  
  output:
  tuple val(meta), path("${meta.integration}.npy")
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import numpy as np

  adata = ad.read_h5ad("${input}")

  ar = adata.obsm["X_emb"]

  np.save("${meta.integration}.npy", ar)
  """

}