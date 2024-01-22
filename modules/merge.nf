process MERGE {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.0"

  publishDir "${params.outdir}", mode: "${params.publish_mode}"

  label "process_high_memory"

  input:
    tuple val(meta), path(original_adata, stageAs: 'input.h5ad')
    tuple val(meta2), path(counts)
    path(obsm)
    path(obs)
  
  output:
    file "merged.h5ad"
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import numpy as np
  import pandas as pd
  import os
  from scipy.sparse import csc_matrix

  adata = ad.read_h5ad("${original_adata}")

  counts_adata = ad.read_h5ad("$counts")
  for layer in counts_adata.layers.keys():
    adata.layers[layer] = csc_matrix(counts_adata.layers[layer]).astype(np.float32)
  adata.X = csc_matrix(counts_adata.X).astype(np.float32)
  del counts_adata

  obsm_paths = "${obsm}".split(' ')

  for obsm_path in obsm_paths:
    array = np.load(obsm_path)
    name = os.path.basename(obsm_path).split('.')[0]
    adata.obsm[name] = np.float32(array)
  
  obs_paths = "${obs}".split(' ')

  for obs_path in obs_paths:
    df = pd.read_pickle(obs_path)
    adata.obs = pd.concat([adata.obs, df], axis=1)
  
  for col in adata.obs.columns:
    if adata.obs[col].dtype == np.float64:
      adata.obs[col] = adata.obs[col].astype(np.float32)

  adata.write('merged.h5ad')
  """
}