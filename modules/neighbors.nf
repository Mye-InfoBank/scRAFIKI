process NEIGHBORS {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*.h5ad")

    script:
    """
    #!/opt/conda/bin/python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.pp.neighbors(adata, use_rep="X_emb")
    adata.write_h5ad("${meta.id}.neighbors.h5ad")
    """
}