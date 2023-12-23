process MERGE_DATASETS {
  container = "bigdatainbiomedicine/sc-python"
  cpus = 4
  memory = {50.GB * task.attempt}
  maxRetries = 4
  errorStrategy = 'retry'

  input:
  path(adatas)
  
  output:
  path("merged.h5ad")
  
  script:
  """
  merge_datasets.py --input ${adatas} --output merged.h5ad
  """
}