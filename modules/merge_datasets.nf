process MERGE_DATASETS {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  path(adatas)
  
  output:
  path("datasets.inner.h5ad"), emit: inner
  path("datasets.outer.h5ad"), emit: outer
  path("batches.txt"), emit: batches

  when:
  task.ext.when == null || task.ext.when
  
  script:
  """
  merge_datasets.py --input ${adatas} --output_inner datasets.inner.h5ad --output_outer datasets.outer.h5ad
  """
}