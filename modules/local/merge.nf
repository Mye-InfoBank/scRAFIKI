process MERGE {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-python"

  publishDir "${params.outdir}", mode: 'symlink'

  label "process_high_memory"

  input:
    tuple val(meta), path(original_adata, stageAs: 'input.h5ad')
    val (integration_names)
    val(integration_embeddings)
    path(integration_adatas)
    tuple val(meta2), path(solo)
    tuple val(meta3), path(decontX)
    val(resolutions)
  
  output:
    file "merged.h5ad"
  
  script:
  """
  #!/usr/bin/env python3

  import anndata as ad
  import pandas as pd

  adata = ad.read_h5ad("${original_adata}")

  for integration_name, integration_embedding, integration_adata_path in zip(
        ["${integration_names.join("\",\"")}"],
        ["${integration_embeddings.join("\",\"")}"],
        ["${integration_adatas.join("\",\"")}"]):
      integration_adata = ad.read_h5ad(integration_adata_path)
      adata.obsm['emb_' + integration_name] = integration_adata.obsm[integration_embedding].copy()
      adata.obsm['X_' + integration_name] = integration_adata.obsm['X_umap'].copy()

      for resolution in ${resolutions}:
        res = float(resolution)
        leiden_key = f"leiden_{res:.2f}"
        adata.obs[integration_name + '_' + leiden_key] = integration_adata.obs[leiden_key].copy()

      del integration_adata

  solo_df = pd.read_csv("$solo", sep="\\t", index_col=0)
  decontX_adata = ad.read_h5ad("$decontX")

  adata.layers["decontX"] = decontX_adata.X
  adata.obs["solo_doublet_score"] = solo_df["doublet"].values
  adata.obs["solo_singlet_score"] = solo_df["singlet"].values
  adata.obs["solo_label"] = solo_df["label"].values

  adata.write('merged.h5ad')
  """
}