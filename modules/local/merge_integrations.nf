process MERGE_INTEGRATIONS {
  container "bigdatainbiomedicine/sc-python"

  label "process_high_memory"

  input:
    tuple val(meta), path(original_adata)
    val  integration_names
    val  integration_types
    file integration_adatas
  
  output:
    file "integrations.h5ad"
  
  script:
  """
  #!/usr/bin/env python3
  import scanpy as sc
  import anndata as ad

  adata = sc.read_h5ad("${original_adata}")

  # Valid integration types are:
  #   - full
  #   - embed
  #   - knn

  for integration_name, integration_type, integration_adata_path in zip(
        ["${integration_names.join("\",\"")}"],
        ["${integration_types.join("\",\"")}"],
        ["${integration_adatas.join("\",\"")}"]):
      integration_adata = ad.read_h5ad(integration_adata_path)

      print(integration_name)
      print(integration_adata)
      
      integration_string = 'X_' + integration_name

      if integration_type == 'full':
        adata.layers[integration_string] = integration_adata.X.copy() # TODO: Check if this works
      elif integration_type == 'embed' or integration_type == 'knn':
        adata.obsm[integration_string] = integration_adata.obsm["X_emb"].copy()
      # elif integration_type == 'knn':
        # TODO: Revisit
        # column_name = 'X_' + integration_name + '_knn'
        # adata.obsp[column_name] = integration_adata.obsp[integration_accession].copy()

      del integration_adata

      # print("Merged integration: " + integration_name)
      # print("Object now looks like:")
      # print(adata)
      # print("")

  adata.write('integrations.h5ad')
  """
}