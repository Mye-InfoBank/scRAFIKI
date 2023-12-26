process INTEGRATE {
  tag "${method}"
  container "bigdatainbiomedicine/sc-python"

  label "process_high"
  label "scale_resources"

  input:
  path(input)
  val(method)
  
  output:
  tuple val(meta), path("${method}.h5ad")
  
  script:
  meta = ["id": "${method}", "integration": "${method}"]
  """
  integrate.py --input ${input} --method ${method} --output ${method}.h5ad
  """
}