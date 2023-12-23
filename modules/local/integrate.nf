process INTEGRATE {
  tag "${method}"
  container "bigdatainbiomedicine/sc-python"

  label "process_medium"
  label "error_retry"

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