process SC_HPL_LEARN {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-schpl:1.0.4"

  publishDir "${params.outdir}/tree", mode: "${params.publish_mode}"
  label "process_high"

  input:
  tuple val(meta), path(adata), val(resolutions), path(tables)
  tuple val(meta2), path(base_tree)

  output:
  tuple val(meta), path("${meta.id}.tree.pkl"), emit: tree
  
  script:
  gpu = task.ext.use_gpu ? '0' : 'None'
  compress = task.ext.compress ? 'True' : 'False'
  """
  #!/usr/bin/env python

  import scHPL
  import anndata as ad
  import pickle
  import pandas as pd
  import numpy as np

  adata = ad.read_h5ad("${adata}")
  adata_latent = ad.AnnData(X=adata.obsm['X_emb'], obs=adata.obs)

  tables = [pd.read_pickle(table) for table in '${tables}'.split(' ')]
  df = pd.concat(tables, axis=1)

  # Prefix all values with their column key
  df = df.apply(lambda x: x.index + '_' + x.astype(str), axis=1)

  resolutions = df.columns

  # For each row, choose a random column
  adata_latent.obs['resolution'] = np.random.choice(resolutions, size=adata_latent.shape[0])
  adata_latent.obs['cluster'] = df.lookup(adata_latent.obs.index, adata_latent.obs['resolution'])

  kwargs = {
      'data': adata_latent,
      'batch_key': 'resolution',
      ${base_tree ? "'tree': pickle.load(open('" + base_tree + "', 'rb'))," : ""}
      'cell_type_key': 'cluster',
      'classifier': 'knn',
      'dynamic_neighbors': True,
      'dimred': False,
      'gpu': ${gpu},
      'compress': ${compress}
  }

  kwargs['batch_order'] = adata_latent.obs['resolution'].unique()

  tree_ref, mp_ref = scHPL.learn_tree(**kwargs)
  pickle.dump(tree_ref, open("${meta.id}.tree.pkl", "wb"))
  """
}

process SC_HPL_PREDICT {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-schpl:1.0.4"

  label "process_high"

  input:
  tuple val(meta), path(adata), path(tree)

  output:
  tuple val(meta), path("${meta.id}.scHPL.pkl")
  
  script:
  gpu = task.ext.use_gpu ? '0' : 'None'
  """
  #!/usr/bin/env python

  import scHPL
  import anndata as ad
  import pickle

  adata = ad.read_h5ad("${adata}")
  tree = pickle.load(open("${tree}", "rb"))

  kwargs = {
      'tree': tree,
      'gpu': ${gpu}
  }

  labels, probabilities = scHPL.predict_labels(adata.obsm['X_emb'], **kwargs)
  adata.obs['${meta.id}_scHPL'] = labels

  df = adata.obs[['${meta.id}_scHPL']]
  df.to_pickle("${meta.id}.scHPL.pkl")
  """
}