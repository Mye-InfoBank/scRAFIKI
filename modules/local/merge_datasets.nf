process MERGE_DATASETS {
  container "bigdatainbiomedicine/sc-python"

  label "process_medium"
  label "error_retry"

  input:
  path(adatas)
  
  output:
  path("merged.h5ad")
  
  script:
  """
  merge_datasets.py --input ${adatas} --output merged.h5ad
  """
}