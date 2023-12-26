process INTEGRATE {
  tag "${method}"
  container "bigdatainbiomedicine/sc-python"

  label "process_high"
  label "scale_resources"

  input:
  tuple val(meta), path(input)
  val(method)
  
  output:
  tuple val(meta_out), path("${method}.h5ad")
  
  script:
  meta_out = ["id": "${method}", "integration": "${method}"]
  """
  integrate.py --input ${input} --method ${method} --output ${method}.h5ad
  """
}