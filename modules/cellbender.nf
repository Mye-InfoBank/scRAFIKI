process CELLBENDER {
  tag "$meta.id"

  container "us.gcr.io/broad-dsde-methods/cellbender:0.3.0"

  label "process_medium"
  label "process_high_memory"

  input:
    tuple val(meta), file(adata)
  
  output:
    tuple val(meta), file("${meta.id}.ambient.h5ad")
  
  script:
  """
  cellbender remove-background \
     --cuda \
     --input ${adata} \
     --output ${meta.id}.ambient.h5ad  
  """
}