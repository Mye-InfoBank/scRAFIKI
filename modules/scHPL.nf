process SC_HPL {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-scarches:1.0"

  label "process_high"

  input:
  tuple val(meta), path(adata)
  tuple val(meta2), path(base_tree)

  output:
  tuple val(meta), path("${meta.id}.tree.pkl"), emit: tree
  
  script:
  tree = base_tree ? "tree=pickle.load(open('${base_tree}', 'rb'))," : ""
  """
  #!/usr/bin/env python

  import scarches as sca
  import scanpy as sc
  import pickle

  adata_full = sc.read_h5ad("${adata}")
  adata_latent = sc.AnnData(X=adata_full.obsm['X_emb'], obs=adata_full.obs)

  adata_latent.obs['celltype_batch'] = adata_latent.obs['cell_type'].astype(str) + '_' + adata_latent.obs['batch'].astype(str)

  batch_counts = adata_latent.obs['batch'].value_counts().to_dict()
  batches_sorted_desc = sorted(batch_counts, key=batch_counts.get, reverse=True)

  tree_ref, mp_ref = sca.classifiers.scHPL.learn_tree(data=adata_latent,
                                                      batch_key='batch',
                                                      batch_order=batches_sorted_desc,
                                                      cell_type_key='celltype_batch',
                                                      classifier='knn',
                                                      ${tree}
                                                      dynamic_neighbors=True,
                                                      dimred=False)
  pickle.dump(tree_ref, open("${meta.id}.tree.pkl", "wb"))
  """
}