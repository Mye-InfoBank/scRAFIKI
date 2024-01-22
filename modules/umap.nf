process UMAP {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.0"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("X_${meta.integration}.npy")

    script:
    """
    #!/opt/conda/bin/python

    import scanpy as sc
    import numpy as np
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.umap(adata)
    
    np.save("X_${meta.integration}.npy", adata.obsm["X_umap"])
    """
}