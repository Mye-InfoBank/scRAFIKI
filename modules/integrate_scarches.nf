process INTEGRATE_SCARCHES {
  tag "${method}"
  container "bigdatainbiomedicine/sc-scarches:1.0"

  label "process_high"

  input:
  tuple val(meta), path(query)
  tuple val(meta2), path(reference), path(scanvi_model, stageAs: "scanvi_model")
  
  output:
  tuple val(meta_out), path("${method}.h5ad"), emit: integrated
  tuple val(meta_out), path("${model}"), emit: model
  
  script:
  method = "scarches"
  model = "model"
  meta_out = ["id": "${method}", "integration": "${method}"]
  """
  #!/usr/bin/env python

  import scarches as sca
  import anndata as ad

  reference_model_path = "${scanvi_model}"
  surgery_model_path = "${model}"
  reference_path = "${reference}"
  query_path = "${query}"

  adata_reference = ad.read_h5ad(reference_path)
  adata_query = ad.read_h5ad(query_path)

  sca.models.SCANVI.prepare_query_anndata(
      adata=adata_query, reference_model=reference_model_path
  )

  surgery_model = sca.models.SCANVI.load_query_data(
      adata_query,
      reference_model_path,
      freeze_dropout=True,
  )

  surgery_epochs = 500
  early_stopping_kwargs_surgery = {
      "early_stopping_monitor": "elbo_train",
      "early_stopping_patience": 10,
      "early_stopping_min_delta": 0.001,
      "plan_kwargs": {"weight_decay": 0.0},
  }

  surgery_model.train(max_epochs=surgery_epochs, **early_stopping_kwargs_surgery)
  surgery_model.save(surgery_model_path, overwrite=True)

  adata_query.obsm["X_emb"] = surgery_model.get_latent_representation(adata_query)
  adata_query.write_h5ad("${method}.h5ad")
  """
}