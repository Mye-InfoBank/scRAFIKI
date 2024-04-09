process SC_HPL {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-scarches:1.0"

  publishDir "${params.outdir}/tree", mode: "${params.publish_mode}"
  label "process_high"

  input:
  tuple val(meta), path(adata), path(tables)

  output:
  tuple val(meta), path("${meta.id}.tree.pkl"), emit: tree
  
  script:
  """
  #!/usr/bin/env python

  import scarches as sca
  import scanpy as sc
  import pickle
  import pandas as pd

  adata = sc.read_h5ad("${adata}")
  adata_latent = sc.AnnData(X=adata.obsm['X_emb'], obs=adata.obs)

  tables = [pd.read_pickle(table) for table in '${tables}'.split(' ')]
  df = pd.concat(tables, axis=1)

  # Prefix all values with their column key
  df = df.apply(lambda x: x.index + '_' + x.astype(str), axis=1)

  # For each row, choose a random column value
  adata_latent.obs['random_cluster'] = df.apply(lambda x: x.sample().values[0], axis=1)
  adata_latent.obs['batch_cluster'] = adata_latent.obs['batch'].astype(str) + '_' + adata_latent.obs['random_cluster'].astype(str)

  kwargs = {
      'data': adata_latent,
      'batch_key': 'batch',
      'cell_type_key': 'batch_cluster',
      'classifier': 'knn',
      'dynamic_neighbors': True,
      'dimred': False
  }

  batch_counts = adata_latent.obs['batch'].value_counts().to_dict()
  batches_sorted_desc = sorted(batch_counts, key=batch_counts.get, reverse=True)
  kwargs['batch_order'] = batches_sorted_desc

  tree_ref, mp_ref = sca.classifiers.scHPL.learn_tree(**kwargs)
  pickle.dump(tree_ref, open("${meta.id}.tree.pkl", "wb"))
  """
}