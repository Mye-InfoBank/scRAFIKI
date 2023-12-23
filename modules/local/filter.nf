process FILTER {
  container = "bigdatainbiomedicine/sc-python"
  cpus = 4
  memory = {50.GB * task.attempt}
  maxRetries = 4
  errorStrategy = 'retry'
  tag "${meta.id}"

  input:
  tuple val(meta), path(input)
  
  output:
  tuple val(meta), path("${meta.id}.filtered.adata")
  
  script:
  min_counts = meta.min_counts ? "--min_counts ${meta.min_counts}" : ""
  max_counts = meta.max_counts ? "--max_counts ${meta.max_counts}" : ""
  min_genes = meta.min_genes ? "--min_genes ${meta.min_genes}" : ""
  max_genes = meta.max_genes ? "--max_genes ${meta.max_genes}" : ""
  max_pct_mito = meta.max_pct_mito ? "--max_pct_mito ${meta.max_pct_mito}" : ""
  """
  filter.py --input ${input} --id ${meta.id} ${min_counts} ${max_counts} ${min_genes} ${max_genes} ${max_pct_mito} --output ${meta.id}.filtered.adata
  """
}