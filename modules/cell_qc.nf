process CELL_QC {
  tag "$meta.id"
  container "bigdatainbiomedicine/sc-rpy:1.2"

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
  n_genes = adata.X.shape[1]
  top_genes = [n for n in [50, 100, 200, 500] if n <= n_genes]
  sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True, percent_top=top_genes)

  qc_columns = ['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts',
        'log1p_total_counts', 'total_counts_mito',
        'log1p_total_counts_mito', 'pct_counts_mito'] + [f"pct_counts_in_top_{n}_genes" for n in top_genes]

  qc = adata.obs[qc_columns].rename(columns={column: f"qc:{column}" for column in qc_columns})

  qc.to_pickle("qc.pkl")
  """
}