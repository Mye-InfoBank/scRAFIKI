process MERGE_INTEGRATIONS {
  input:
    file original_adata
    file integration_adatas
    value integration_names
  
  output:
    file "integrations.h5ad"
  
  script:
  """
  #!/usr/bin/env python
  import scanpy as sc

  adata = sc.read_h5ad(original_adata)

  for integration_name, integration_adata in zip(integration_names, integration_adatas):
      integration_adata = sc.read_h5ad(integration_adata)
      column_name = 'X_' + integration_name
      adata.obsm[column_name] = integration_adata.obsm[column_name]

  adata.write('integrations.h5ad')
  """
}