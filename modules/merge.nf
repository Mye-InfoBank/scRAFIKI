process MERGE {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.0"

  publishDir "${params.outdir}", mode: "${params.publish_mode}"

  label "process_high_memory"

  input:
    tuple val(meta), path(base, stageAs: 'base.h5ad')
    path(obsm)
    path(obs)
  
  output:
    file "merged.h5ad"
    file "metadata.pkl"
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import numpy as np
  import pandas as pd
  import os
  from scipy.sparse import csc_matrix

  adata = ad.read_h5ad("${base}")
  adata.obsm = {}

  for layer in adata.layers.keys():
    adata.layers[layer] = csc_matrix(adata.layers[layer]).astype(np.float32)
  adata.X = csc_matrix(adata.X).astype(np.float32)

  obsm_paths = "${obsm}".split(' ')

  for obsm_path in obsm_paths:
    df = pd.read_pickle(obsm_path)
    df = df.reindex(adata.obs_names)
    name = os.path.basename(obsm_path).split('.')[0]
    adata.obsm[name] = np.float32(df.to_numpy())
  
  obs_paths = "${obs}".split(' ')

  for obs_path in obs_paths:
    df = pd.read_pickle(obs_path)
    df = df.reindex(adata.obs_names)
    adata.obs = pd.concat([adata.obs, df], axis=1)

  for col in adata.obs.columns:
    if adata.obs[col].dtype == np.float64:
      adata.obs[col] = adata.obs[col].astype(np.float32)

  integration_order = ["scarches", "scanvi", "scvi", "scgen", "desc",
                       "bbknn", "combat", "harmony", "mnn",
                       "scanorama", "trvaep", "unintegrated"]

  current_index = 1
  for integration_method in integration_order:
    integration_key = f"X_{integration_method}"
    if integration_key in adata.obsm.keys():
      target_key = f"X_{current_index}_{integration_method}"
      adata.obsm[target_key] = adata.obsm[integration_key]
      del adata.obsm[integration_key]
      current_index += 1

  adata.write('merged.h5ad')
  adata.obs.to_pickle('metadata.pkl')
  """
}