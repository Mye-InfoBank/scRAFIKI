process INTEGRATE {
  container = "bigdatainbiomedicine/sc-python"
  cpus = 4
  memory = {50.GB * task.attempt}
  maxRetries = 4
  errorStrategy = 'retry'
  tag "${method}"

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