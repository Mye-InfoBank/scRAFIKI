process MERGE_DATASETS {
  container "bigdatainbiomedicine/sc-rpy"

  label "process_medium"

  input:
  path(adatas)
  
  output:
  path("merged_datasets.h5ad")
  
  script:
  """
  merge_datasets.py --input ${adatas} --output merged_datasets.h5ad
  """
}