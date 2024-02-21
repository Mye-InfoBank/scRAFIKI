process EXTRACT_EMBEDDING {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  tuple val(meta), path(input)
  val(embedding)
  
  output:
  tuple val(meta), path("${meta.integration}.pkl")
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import pandas as pd

  adata = ad.read_h5ad("${input}")

  ar = adata.obsm["${embedding}"]

  df = pd.DataFrame(ar, index=adata.obs_names)
  df.to_pickle("${meta.integration}.pkl")
  """

}