process UMAP {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.0"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("X_${meta.id}.pkl")

    script:
    """
    #!/opt/conda/bin/python

    import scanpy as sc
    import pandas as pd
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.umap(adata)
    
    df = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names)
    df.to_pickle("X_${meta.id}.pkl")
    """
}