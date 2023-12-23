process INTEGRATE {
  tag "${method}"
  container "bigdatainbiomedicine/sc-python"

  label "process_high"
  label "scale_resources"

  input:
  path(input)
  val(method)
  
  output:
  path("${method}.adata")
  
  script:
  """
  integrate.py --input ${input} --method ${method} --output ${method}.adata
  """
}