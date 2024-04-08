process SC_HPL {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-scarches:1.0"

  label "process_high"

  input:
  tuple val(meta ), path(adata)
  tuple val(meta2), path(base_adata)
  tuple val(meta3), path(base_tree)

  output:
  tuple val(meta), path("${meta.id}.tree.pkl"), emit: tree
  
  script:
  """
  #!/usr/bin/env python

  import scarches as sca
  import scanpy as sc
  import pickle

  has_base = ${base_adata && base_tree ? "True" : "False"}

  adata_full = sc.read_h5ad("${adata}")
  adata_latent = sc.AnnData(X=adata_full.obsm['X_emb'], obs=adata_full.obs)

  adata_latent.obs['celltype_batch'] = adata_latent.obs['cell_type'].astype(str) + '_' + adata_latent.obs['batch'].astype(str)

  kwargs = {
      'data': adata_latent,
      'batch_key': 'batch',
      'cell_type_key': 'celltype_batch',
      'classifier': 'knn',
      'dynamic_neighbors': True,
      'dimred': False
  }

  batch_counts = adata_latent.obs['batch'].value_counts().to_dict()
  batches_sorted_desc = sorted(batch_counts, key=batch_counts.get, reverse=True)

  if has_base:
    adata_base = sc.read_h5ad("${base_adata}")
    batches_base = adata_base.obs['batch'].unique().tolist()
    batches_sorted_desc = [x for x in batches_sorted_desc if x not in batches_base]

    tree_ref = pickle.load(open("${base_tree}", "rb"))
    kwargs['tree'] = tree_ref
  
  kwargs['batch_order'] = batches_sorted_desc

  tree_ref, mp_ref = sca.classifiers.scHPL.learn_tree(**kwargs)
  pickle.dump(tree_ref, open("${meta.id}.tree.pkl", "wb"))
  """
}