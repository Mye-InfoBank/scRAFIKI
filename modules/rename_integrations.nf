process RENAME_INTEGRATIONS {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.0"

  publishDir "${params.outdir}", mode: "${params.publish_mode}"

  label "process_high_memory"

  input:
    tuple val(meta), path(adata)
  
  output:
    file "merged.h5ad"
    file "metadata.pkl"
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad

  adata = ad.read_h5ad("${adata}")

  integration_order = ["scarches", "scanvi", "scvi", "scgen", "desc",
                       "bbknn", "combat", "harmony", "mnn",
                       "scanorama", "trvaep", "unintegrated"]

  name_conversion = {}

  current_index = 1
  for integration_method in integration_order:
    integration_key = f"X_{integration_method}"
    if integration_key in adata.obsm.keys():
      target_key = f"X_{current_index}_{integration_method}"
      adata.obsm[target_key] = adata.obsm[integration_key]
      del adata.obsm[integration_key]
      name_conversion[integration_key] = target_key
      current_index += 1

  adata.uns['integration_names'] = name_conversion

  adata.write('atlas.h5ad')
  """
}