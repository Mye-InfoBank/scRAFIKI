process MERGE_DATASETS {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  path(adatas)
  
  output:
  path("merged_datasets.h5ad")

  when:
  task.ext.when == null || task.ext.when
  
  script:
  """
  merge_datasets.py --input ${adatas} --output merged_datasets.h5ad
  """
}