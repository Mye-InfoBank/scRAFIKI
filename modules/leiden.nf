process LEIDEN {
    tag "${meta.id}:${resolution}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_medium"

    input:
    tuple val(meta), path(adata), val(resolution)

    output:
    tuple val(meta), val("${clustering_key}"), path("*.h5ad")

    script:
    clustering_key = "leiden_${resolution}"
    """
    #!/opt/conda/bin/python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.leiden(adata, resolution=${resolution}, key_added="${clustering_key}")
    adata.write_h5ad("${meta.id}.${clustering_key}.h5ad")
    """
}