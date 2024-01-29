process SOLO {
    tag "${new_meta.id}"
    container "bigdatainbiomedicine/sc-scib:1.0"

    label "process_medium"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(scvi_model)
        each batch

    output:
        tuple val(new_meta), path("${new_meta.id}.solo.pkl")

    script:
    new_meta = meta + [id: "${meta.id}-${batch}", batch: batch]
    """
    #!/usr/bin/env python3

    import scvi
    import scanpy as sc
    from threadpoolctl import threadpool_limits
    threadpool_limits(${task.cpus})

    adata = sc.read_h5ad("${adata}")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    scvi.model.SCANVI.setup_anndata(adata_hvg, batch_key="batch", labels_key="cell_type", unlabeled_category="Unknown")
    scvi_model = scvi.model.SCANVI.load("${scvi_model}", adata=adata_hvg)
    solo = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch="${batch}")
    solo.train()
    res = solo.predict()
    res["doublet_label"] = solo.predict(False)

    res.to_pickle("${new_meta.id}.solo.pkl")
    """
}

process MERGE_SOLO {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-scib:1.0"

    label "process_medium"

    input:
        tuple val(meta), path(adata)
        path(solo_results)

    output:
        tuple val(meta), path("solo.pkl")

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import anndata as ad

    adata = ad.read_h5ad("${adata}")

    solo_paths = "${solo_results}".split(" ")
    solo_results = [pd.read_pickle(path) for path in solo_paths]

    solo_results = pd.concat(solo_results)

    # Reorder the cells to match the original adata
    solo_results = solo_results.reindex(adata.obs_names)

    solo_results.to_pickle("solo.pkl")
    """
}