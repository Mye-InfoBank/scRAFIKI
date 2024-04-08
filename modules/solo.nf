process SOLO {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-scib:1.0"

    label "process_medium"
    label "process_high_memory"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(scvi_model)
        val(has_celltypes)
        val(batch)

    output:
        tuple val(new_meta), path("${new_meta.id}.solo.pkl")

    script:
    new_meta = [id: "${batch}"]
    """
    #!/usr/bin/env python3

    import scvi
    import scanpy as sc
    import pandas as pd
    from threadpoolctl import threadpool_limits
    threadpool_limits(${task.cpus})

    adata = sc.read_h5ad("${adata}")
    adata_batch = adata[adata.obs.batch == "${batch}"]

    model_path = "${scvi_model}"

    if ${has_celltypes ? "True" : "False"}:
        # scvi.model.SCANVI.setup_anndata(adata, batch_key="batch", labels_key="cell_type", unlabeled_category="Unknown")
        scvi.model.SCANVI.prepare_query_anndata(
            adata=adata, reference_model=model_path
        )
        scvi_model = scvi.model.SCANVI.load(model_path, adata=adata)
    else:
        # scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
        scvi.model.SCVI.prepare_query_anndata(
            adata=adata, reference_model=model_path
        )
        scvi_model = scvi.model.SCVI.load(model_path, adata=adata)

    solo = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch="${batch}")

    minibatch_size = min(128, len(adata_batch))
    worked = False
    while not worked and minibatch_size > 100:
        try:
            solo.train(batch_size=minibatch_size)
            worked = True
        except ValueError:
            print("Minibatch size did not work, trying again with smaller minibatch size")
            minibatch_size -= 1

    if worked:
        solo_res = solo.predict()
        solo_res["doublet_label"] = solo.predict(False)
    else:
        solo_res = pd.DataFrame(index=adata_batch.obs.index)
        solo_res["doublet_label"] = "Unknown"

    solo_res.to_pickle("${new_meta.id}.solo.pkl")
    """
}