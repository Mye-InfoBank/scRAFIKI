process NEIGHBORS {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*.h5ad"), emit: adata

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

process UMAP {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*.h5ad"), emit: adata

    script:
    """
    #!/opt/conda/bin/python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.umap(adata)
    adata.write_h5ad("${meta.id}.umap.h5ad")
    """
}

process LEIDEN {
    tag "${meta.id}:${resolution}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_single"

    input:
    tuple val(meta), path(adata)
    each resolution

    output:
    tuple val(meta), val(resolution), path("*.h5ad"), emit: adata

    script:
    """
    #!/opt/conda/bin/python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.leiden(adata, resolution=${resolution})
    adata.write_h5ad("${meta.id}.res_${resolution}.leiden.h5ad")
    """
}

process CELLTYPIST_MAJORITY {
  tag "${meta.id}"

  container "bigdatainbiomedicine/sc-rpy"
  label "process_medium"

  input:
  tuple val(meta), val(resolution), path(adata)
  tuple val(meta2), path(celltypist)
  
  output:
  tuple val(meta), val(resolution), path("${meta.id}.res_${resolution}.majority.h5ad"), emit: adata
  
  script:
  """
  celltypist_majority.py --input_clustering ${adata} --input_celltypist ${celltypist} --output ${meta.id}.res_${resolution}.majority.h5ad
  """
}

process MERGE_UMAP_LEIDEN {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_medium"


    input:
    tuple val(meta), path(adata_umap), val(leiden_resolutions), path(adata_leiden)

    output:
    tuple val(meta), path("*.h5ad"), emit: adata

    script:
    """
    #!/opt/conda/bin/python

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
        adata_umap.obs[f"leiden_{res:.2f}_celltypist_majority"] = tmp_adata.obs["celltypist_majority"]
    adata_umap.write_h5ad("${meta.id}.umap_leiden.h5ad")
    """
}



workflow NEIGHBORS_LEIDEN_UMAP {
    take:
    adata
    leiden_res
    ch_celltypist

    main:
    NEIGHBORS(adata)
    UMAP(NEIGHBORS.out.adata)
    LEIDEN(NEIGHBORS.out.adata, leiden_res)
    CELLTYPIST_MAJORITY(LEIDEN.out.adata, ch_celltypist)

    MERGE_UMAP_LEIDEN(
        UMAP.out.adata.join(
            CELLTYPIST_MAJORITY.out.adata.groupTuple(by: 0), by: 0
        )
    )

    emit:
    adata = MERGE_UMAP_LEIDEN.out.adata
}
