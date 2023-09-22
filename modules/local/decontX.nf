process DECONTX {
  input:
    tuple val(batch), file(adata)
  
  output:
    tuple val(batch), file("*_ambient.h5ad")
  
  script:
  """
  decontX.R ${adata}
  """
}