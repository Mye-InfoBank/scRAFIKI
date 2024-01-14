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

process UMAP {
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
    sc.tl.umap(adata)
    adata.write_h5ad("${meta.id}.umap.h5ad")
    """
}

process LEIDEN {
    tag "${meta.id}:${resolution}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_medium"

    input:
    tuple val(meta), path(adata)
    each resolution

    output:
    tuple val(meta), val(resolution), path("*.h5ad")

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

process ENTROPY {
  tag "${meta.id}"

  container "bigdatainbiomedicine/sc-rpy"
  label "process_medium"

  input:
  tuple val(meta), val(resolution), path(adata)
  
  output:
  tuple val(meta), val(resolution), path("${meta.id}.res_${resolution}.entropy.h5ad")
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import anndata2ri
  import rpy2.robjects as ro
  rogue = ro.packages.importr('ROGUE')
  _ = ro.packages.importr('tibble')

  adata = ad.read_h5ad("$adata")
  sce = anndata2ri.py2rpy(adata)
  expression = ro.r("counts")(sce)

  res = float("${resolution}")

  label_col = f"leiden_{res:.2f}"
  sample_col = "batch"

  labels = anndata2ri.py2rpy(adata.obs[label_col])
  samples = anndata2ri.py2rpy(adata.obs[sample_col])

  result = rogue.rogue(expression, labels=labels, samples=samples, platform="UMI")
  result = anndata2ri.rpy2py(result)
  entropy_dict = result.to_dict()

  adata.obs[f"{label_col}_entropy"] = adata.obs.apply(
    lambda row: entropy_dict[row[label_col]][row[sample_col]], axis=1
  )

  adata.write_h5ad("${meta.id}.res_${resolution}.entropy.h5ad")
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
  tuple val(meta), val(resolution), path("${meta.id}.res_${resolution}.majority.h5ad")
  
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
    tuple val(meta), path("*.h5ad")

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



workflow CLUSTERING {
    take:
        ch_adata
        ch_resolutions
        ch_celltypist

    main:
        NEIGHBORS(ch_adata)
        UMAP(NEIGHBORS.out)
        LEIDEN(NEIGHBORS.out, ch_resolutions)
        ENTROPY(LEIDEN.out)
        CELLTYPIST_MAJORITY(LEIDEN.out, ch_celltypist)

        MERGE_UMAP_LEIDEN(
            UMAP.out.join(
                CELLTYPIST_MAJORITY.out.groupTuple(by: 0), by: 0
            )
        )

    emit:
        MERGE_UMAP_LEIDEN.out
}
