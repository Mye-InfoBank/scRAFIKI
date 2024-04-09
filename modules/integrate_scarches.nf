process INTEGRATE_SCARCHES {
  tag "${method}"
  container "bigdatainbiomedicine/sc-scib:1.3"

  label "process_high"

  input:
  tuple val(meta), path(query)
  tuple val(meta2), path(base_model, stageAs: "base_model")
  val(has_celltypes)
  
  output:
  tuple val(meta_out), path("${method}.integrated.h5ad"), emit: integrated
  tuple val(meta_out), path("${model}"), emit: model
  
  script:
  method = "scarches"
  model = "model"
  meta_out = ["id": "${method}", "integration": "${method}"]
  """
  #!/usr/bin/env python

  import scarches as sca
  import anndata as ad
  import pandas as pd

  reference_model_path = "${base_model}"
  surgery_model_path = "${model}"
  query_path = "${query}"

  adata_query = ad.read_h5ad(query_path)

  adata_output = adata_query.copy()

  if ${has_celltypes ? "True" : "False"}:
    sca.models.SCANVI.prepare_query_anndata(
        adata=adata_query, reference_model=reference_model_path
    )

    surgery_model = sca.models.SCANVI.load_query_data(
        adata_query,
        reference_model_path,
        freeze_dropout=True,
    )
  else:
    sca.models.SCVI.prepare_query_anndata(
        adata=adata_query, reference_model=reference_model_path
    )

    surgery_model = sca.models.SCVI.load_query_data(
        adata_query,
        reference_model_path,
        freeze_dropout=True,
    )

  surgery_epochs = 100
  early_stopping_kwargs_surgery = {
      "early_stopping_monitor": "elbo_train",
      "early_stopping_patience": 10,
      "early_stopping_min_delta": 0.001,
      "plan_kwargs": {"weight_decay": 0.0},
  }

  surgery_model.train(max_epochs=surgery_epochs, **early_stopping_kwargs_surgery)
  surgery_model.save(surgery_model_path, overwrite=True)

  adata_output.obsm["X_emb"] = surgery_model.get_latent_representation()
  adata_output.write_h5ad("${method}.integrated.h5ad")
  """
}