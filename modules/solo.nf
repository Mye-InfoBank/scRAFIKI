process SOLO {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-scib:1.0"

    label "process_medium"
    label "process_high_memory"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(scvi_model)
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

    scvi.model.SCANVI.setup_anndata(adata, batch_key="batch", labels_key="cell_type", unlabeled_category="Unknown")
    scvi_model = scvi.model.SCANVI.load("${scvi_model}", adata=adata)

    solo = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch="${batch}")

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

    solo_res = solo.predict()
    solo_res["doublet_label"] = solo.predict(False)

    solo_res.to_pickle("${new_meta.id}.solo.pkl")
    """
}