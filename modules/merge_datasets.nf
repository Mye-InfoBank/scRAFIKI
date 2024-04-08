process MERGE_DATASETS {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  path(adatas)
  val(min_cells)
  val(custom_genes)
  
  output:
  path("datasets.union.h5ad")       , emit: union
  path("datasets.intersection.h5ad"), emit: intersection

  when:
  task.ext.when == null || task.ext.when
  
  script:
  """
  merge_datasets.py --input ${adatas} --custom_genes ${custom_genes.join(" ")} --min_cells ${min_cells} --output_intersection datasets.intersection.h5ad --output_union datasets.union.h5ad
  """
}