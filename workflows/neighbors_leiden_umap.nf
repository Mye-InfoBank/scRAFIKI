process NEIGHBORS {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-python"
    label "process_medium"

    input:
    tuple val(meta), path(adata)
    each use_rep

    output:
    tuple val(meta), val(use_rep), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.pp.neighbors(adata, use_rep="${use_rep}")
    adata.write_h5ad("${meta.id}.rep_${use_rep}.neighbors.h5ad")
    """
}

process UMAP {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-python"
    label "process_medium"

    input:
    tuple val(meta), val(use_rep), path(adata)

    output:
    tuple val(meta), val(use_rep), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.umap(adata)
    adata.write_h5ad("${meta.id}.rep_${use_rep}.umap.h5ad")
    """
}

process LEIDEN {
    tag "${meta.id}:${resolution}"

    container = "bigdatainbiomedicine/sc-python"
    label "process_single"

    input:
    tuple val(meta), val(use_rep), path(adata)
    each resolution

    output:
    tuple val(meta), val(use_rep), val(resolution), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.leiden(adata, resolution=${resolution})
    adata.write_h5ad("${meta.id}.res_${resolution}.rep_${use_rep}.leiden.h5ad")
    """
}

process MERGE_UMAP_LEIDEN {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-python"
    label "process_medium"


    input:
    tuple val(meta), val(use_rep), path(adata_umap), val(leiden_resolutions), path(adata_leiden)

    output:
    tuple val(meta), val(use_rep), path("*.h5ad"), emit: adata

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
    adata_umap.write_h5ad("${meta.id}.rep_${use_rep}.umap_leiden.h5ad")
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

    MERGE_UMAP_LEIDEN(
        UMAP.out.adata.join(
            LEIDEN.out.groupTuple(by: [0, 1]), by: [0, 1]
        )
    )

    emit:
    adata = MERGE_UMAP_LEIDEN.out.adata
}
