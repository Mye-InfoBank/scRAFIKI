process CELL_QC {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"

  input:
  tuple val(meta), path(adata)
  
  output:
  tuple val(meta), path("qc.pkl")

  when:
  task.ext.when == null || task.ext.when
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import scanpy as sc

  adata = ad.read_h5ad("${adata}")

  adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")
  sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)

  qc_columns = ['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts',
        'log1p_total_counts', 'pct_counts_in_top_50_genes',
        'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes',
        'pct_counts_in_top_500_genes', 'total_counts_mito',
        'log1p_total_counts_mito', 'pct_counts_mito']

  qc = adata.obs[qc_columns].rename(columns={column: f"qc:{column}" for column in qc_columns})

  qc.to_pickle("qc.pkl")
  """
}