process DECONTX {
  input:
    tuple file(adata), val(batch)
  
  output:
    file "${batch}_ambient.h5ad"
  
  script:
  """
  decontX.R ${adata} ${batch}
  """
}