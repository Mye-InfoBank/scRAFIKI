process SC_HPL {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-schpl:1.0.1"

  publishDir "${params.outdir}/tree", mode: "${params.publish_mode}"
  label "process_high"

  input:
  tuple val(meta), path(adata), val(resolutions), path(tables)

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

  adata = ad.read_h5ad("${adata}")
  adata_latent = ad.AnnData(X=adata.obsm['X_emb'], obs=adata.obs)

  tables = zip([${resolutions.join(", ")}], [pd.read_pickle(table) for table in '${tables}'.split(' ')])
  tables = [x[1] for x in sorted(tables, key=lambda x: x[0])]
  df = pd.concat(tables, axis=1)

  adatas = []

  for column in df.columns:
    adata_resolution = adata_latent.copy()
    adata_resolution.obs['cluster'] = column + '_' + df[column].astype(str)
    adatas.append(adata_resolution)

  adata_clusterings = ad.concat(adatas, keys=df.columns, label='resolution', index_unique='-')

  kwargs = {
      'data': adata_clusterings,
      'batch_key': 'resolution',
      'cell_type_key': 'cluster',
      'classifier': 'knn',
      'dynamic_neighbors': True,
      'dimred': False,
      'print_tree': False,
      'gpu': ${gpu},
      'compress': ${compress}
  }

  kwargs['batch_order'] = df.columns

  tree_ref, mp_ref = scHPL.learn_tree(**kwargs)
  pickle.dump(tree_ref, open("${meta.id}.tree.pkl", "wb"))
  """
}