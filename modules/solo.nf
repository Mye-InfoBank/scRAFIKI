process SOLO {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-scib:1.0"

    label "process_medium"
    label "process_high_memory"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(scvi_model)
        val batches

    output:
        tuple val(meta), path("${meta.id}.solo.pkl")

    script:
    """
    #!/usr/bin/env python3

    import scvi
    import scanpy as sc
    import pandas as pd
    from threadpoolctl import threadpool_limits
    threadpool_limits(${task.cpus})

    adata = sc.read_h5ad("${adata}")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    scvi.model.SCANVI.setup_anndata(adata_hvg, batch_key="batch", labels_key="cell_type", unlabeled_category="Unknown")
    scvi_model = scvi.model.SCANVI.load("${scvi_model}", adata=adata_hvg)

    results = []

    batches = "${batches.join(" ")}".split(" ")
    for batch in batches:
        solo = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch=batch)

        minibatch_size = 128
        worked = False
        while not worked and minibatch_size > 100:
            try:
                solo.train(batch_size=minibatch_size)
                worked = True
            except ValueError:
                print("Minibatch size did not work, trying again with smaller minibatch size")
                minibatch_size -= 1
                pass

        batch_res = solo.predict()
        batch_res["doublet_label"] = solo.predict(False)

        results.append(batch_res)
    
    solo_res = pd.concat(results)

    # Reorder the cells to match the original adata
    solo_res = solo_res.reindex(adata.obs_names)

    solo_res.to_pickle("${meta.id}.solo.pkl")
    """
}