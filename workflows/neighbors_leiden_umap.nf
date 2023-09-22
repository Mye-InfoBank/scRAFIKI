process NEIGHBORS {
    container = "bigdatainbiomedicine/sc-python"
    cpus 8
    memory {50.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(adata)
    val use_rep

    output:
    path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.pp.neighbors(adata, use_rep="${use_rep}")
    adata.write_h5ad("neighbors.h5ad")
    """
}

process UMAP {
    container = "bigdatainbiomedicine/sc-python"
    cpus 8
    memory {50.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(adata)

    output:
    path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.umap(adata)
    adata.write_h5ad("umap.h5ad")
    """
}

process LEIDEN {
    container = "bigdatainbiomedicine/sc-python"
    cpus 1
    memory {50.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(adata)
    each resolution

    output:
    tuple val(resolution), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.leiden(adata, resolution=${resolution})
    adata.write_h5ad("res_${resolution}.leiden.h5ad")
    """
}

process MERGE_UMAP_LEIDEN {
    container = "bigdatainbiomedicine/sc-python"
    cpus 1
    memory {50.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(adata_umap), val(leiden_resolutions), path(adata_leiden)

    output:
    path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    resolutions = ${leiden_resolutions}
    if not isinstance(resolutions, list):
        resolutions = [resolutions]
    leiden_adatas = "${adata_leiden}".split(" ")

    adata_umap = sc.read_h5ad("${adata_umap}")
    for res, adata_path in zip(resolutions, leiden_adatas):
        tmp_adata = sc.read_h5ad(adata_path)
        adata_umap.obs[f"leiden_{res:.2f}"] = tmp_adata.obs["leiden"]
    adata_umap.write_h5ad("umap_leiden.h5ad")
    """
}



workflow NEIGHBORS_LEIDEN_UMAP {
    take:
    adata
    neihbors_rep
    leiden_res

    main:
    NEIGHBORS(adata, neihbors_rep)
    UMAP(NEIGHBORS.out.adata)
    LEIDEN(NEIGHBORS.out.adata, leiden_res)

    MERGE_UMAP_LEIDEN(UMAP.out.adata.combine(LEIDEN.out.collect()))

    emit:
    adata = MERGE_UMAP_LEIDEN.out.adata
}
