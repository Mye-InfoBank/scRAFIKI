process CLEAN_ADATA {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_high_memory"

  input:
    tuple val(meta), path(adata)
    val(integration)
  
  output:
    tuple val(meta), path("cleaned.h5ad"), emit: adata
    tuple val(meta), path("embedding.pkl"), emit: embedding
    tuple val(meta), path("X_Global.pkl"), emit: umap
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import pandas as pd

  adata = ad.read_h5ad("${adata}")
  integration = "${integration}"

  embedding = adata.obsm[integration]
  umap_accessor = [obsm for obsm in adata.obsm if obsm.startswith("X_") and integration in obsm][0]
  umap = adata.obsm[umap_accessor]

  adata.obsm = {"X_emb": embedding}

  drop_obs = [obs for obs in adata.obs.columns if "_leiden_" in obs and not obs.startswith(integration)]
  adata.obs = adata.obs.drop(columns=drop_obs)

  rename_obs = [obs for obs in adata.obs.columns if "_leiden_" in obs and obs.startswith(integration)]
  for obs in rename_obs:
    adata.obs.rename(columns={obs: obs.replace(f"{integration}_", "global_")}, inplace=True)

  pd.DataFrame(embedding, index=adata.obs_names).to_pickle("embedding.pkl")
  pd.DataFrame(umap, index=adata.obs_names).to_pickle("X_Global.pkl")

  adata.write("cleaned.h5ad")
  """
}