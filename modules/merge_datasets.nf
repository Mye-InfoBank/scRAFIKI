process MERGE_DATASETS {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  path(adatas)
  tuple val(meta), path(base)
  val(min_cells)
  val(custom_genes)
  
  output:
  path("datasets.union.h5ad")       , emit: union
  path("datasets.intersection.h5ad"), emit: intersection
  path("datasets.transfer.h5ad")    , emit: transfer, optional: true

  when:
  task.ext.when == null || task.ext.when
  
  script:
  if (base) {
    """
    merge_datasets.py --input ${adatas} \\
                      --base ${base} \\
                      --custom_genes ${custom_genes.join(" ")} \\
                      --min_cells ${min_cells} \\
                      --output_intersection datasets.intersection.h5ad \\
                      --output_union datasets.union.h5ad \\
                      --output_transfer datasets.transfer.h5ad
    """
  } else {
    """
    merge_datasets.py --input ${adatas} \\
                      --custom_genes ${custom_genes.join(" ")} \\
                      --min_cells ${min_cells} \\
                      --output_intersection datasets.intersection.h5ad \\
                      --output_union datasets.union.h5ad
    """
  }
}