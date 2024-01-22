#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process SOLO {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-scib:1.0"

    label "process_medium"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(scvi_model)

    output:
        tuple val(meta), path("${meta.id}.solo.pkl")

    script:
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
    solo = scvi.external.SOLO.from_scvi_model(scvi_model)
    solo.train()
    res = solo.predict()
    res["label"] = solo.predict(False)

    res.to_pickle("${meta.id}.solo.pkl")
    """
}