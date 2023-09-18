process DECONTX {
  input:
    file adata
    val batch
  
  output:
    file "${batch}_ambient.h5ad"
  
  script:
  """
  decontX.R ${adata} ${batch}
  """
}