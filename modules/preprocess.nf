process PREPROCESS {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.2"

  label "process_medium"

  input:
  tuple val(meta), path(input)
  val(custom_metadata)
  
  output:
  tuple val(meta), path("${meta.id}.preprocessed.h5ad"), emit: adata, optional: true
  tuple val(meta), path("${meta.id}.problems.txt"), emit: problems, optional: true
  
  script:
  min_counts = meta.min_counts ? "--min_counts ${meta.min_counts}" : ""
  max_counts = meta.max_counts ? "--max_counts ${meta.max_counts}" : ""
  min_genes = meta.min_genes ? "--min_genes ${meta.min_genes}" : ""
  max_genes = meta.max_genes ? "--max_genes ${meta.max_genes}" : ""
  max_pct_mito = meta.max_pct_mito ? "--max_pct_mito ${meta.max_pct_mito}" : ""
  no_symbols = meta.no_symbols && meta.no_symbols.toLowerCase() == "true" ? "--no-symbols" : ""
  sure_raw = meta.sure_raw && meta.sure_raw.toLowerCase() == "true" ? "--sure_raw" : ""
  """
  preprocess.py --input ${input} --custom_metadata ${custom_metadata.join(" ")} --id ${meta.id} ${sure_raw} ${no_symbols} ${min_counts} ${max_counts} ${min_genes} ${max_genes} ${max_pct_mito} --output ${meta.id}.preprocessed.h5ad --problems ${meta.id}.problems.txt
  """
}