process H5_TO_H5AD {
  tag "${meta.id}"

  container "bigdatainbiomedicine/sc-rpy:1.0"
  label "process_medium"

  input:
  tuple val(meta), path(h5)
  
  output:
  tuple val(meta), path("${meta.id}.h5ad")
  
  script:
  """
  #!/opt/conda/bin/python

  import scanpy as sc

  adata = sc.read_10x_h5("${h5}")

  adata.write_h5ad("${meta.id}.h5ad")
  """
}