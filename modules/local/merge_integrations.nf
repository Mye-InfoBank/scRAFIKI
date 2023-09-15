process MERGE_INTEGRATIONS {
  input:
    file original_adata
    file integration_adatas
    val integration_names
  
  output:
    file "integrations.h5ad"
  
  script:
  """
  #!/usr/bin/python3
  import scanpy as sc

  adata = sc.read_h5ad("${original_adata}")

  for integration_name, integration_adata_path in zip(["${integration_names.join("\",\"")}"], ["${integration_adatas.join("\",\"")}"]):
      integration_adata = sc.read_h5ad(integration_adata_path)
      column_name = 'X_' + integration_name
      adata.obsm[column_name] = integration_adata.obsm[column_name]
      del integration_adata

  adata.write('integrations.h5ad')
  """
}