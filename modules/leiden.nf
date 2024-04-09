process LEIDEN {
    tag "${meta.id}:${resolution}"

    container = "bigdatainbiomedicine/sc-rpy:1.2"
    label "process_medium"

    input:
    tuple val(meta), path(adata), val(resolution)

    output:
    tuple val(new_meta), path("${new_meta.id}.clustering.h5ad"), emit: adata
    tuple val(new_meta), path("${new_meta.id}.clustering.pkl"), emit: table

    script:
    new_meta = meta + [id: "${meta.id}_leiden_${resolution}"]
    """
    #!/opt/conda/bin/python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.leiden(adata, resolution=${resolution}, key_added="${new_meta.id}")

    adata.obs["${new_meta.id}"].to_pickle("${new_meta.id}.clustering.pkl")
    adata.write_h5ad("${new_meta.id}.clustering.h5ad")
    """
}