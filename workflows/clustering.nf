process NEIGHBORS {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rapids:0.5"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*.h5ad")

    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import rapids_singlecell as rsc

    adata = ad.read_h5ad("${adata}")
    rsc.utils.anndata_to_GPU(adata)
    rsc.pp.neighbors(adata, use_rep="X_emb")
    adata.write_h5ad("${meta.id}.neighbors.h5ad")
    """
}

process UMAP {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.0"
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

    container = "bigdatainbiomedicine/sc-rpy:1.0"
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
  tag "${meta.id}:${resolution}"

  container "bigdatainbiomedicine/sc-rpy:1.0"
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
  import numpy as np
  rogue = ro.packages.importr('ROGUE')
  _ = ro.packages.importr('tibble')

  adata = ad.read_h5ad("$adata")
  sce = anndata2ri.py2rpy(adata)
  expression = ro.r("assay")(sce, "counts")

  res = float("${resolution}")

  label_col = "leiden"
  sample_col = "batch"

  labels = anndata2ri.py2rpy(adata.obs[label_col].astype(str))
  samples = anndata2ri.py2rpy(adata.obs[sample_col])

  smoothness = 0.5

  while not "entropy" in adata.obs.columns:
    try:
      result = rogue.rogue(expression, labels=labels, samples=samples, platform="UMI", span=smoothness)
      result = anndata2ri.rpy2py(result)
      entropy_dict = result.to_dict()

      adata.obs["entropy"] = adata.obs.apply(
        lambda row: entropy_dict.get(row[label_col], {}).get(row[sample_col], np.nan), axis=1
      )
    except:
      print(f"Failed to compute entropy for ${meta.id} at resolution ${resolution} with smoothness {smoothness}")
      
      if smoothness > 5:
        print("Giving up")
        break
      else:
        smoothness = smoothness * 1.2

  adata.write_h5ad("${meta.id}.res_${resolution}.entropy.h5ad")
  """
}

process CELLTYPIST_MAJORITY {
  tag "${meta.id}"

  container "bigdatainbiomedicine/sc-rpy:1.0"
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

    container = "bigdatainbiomedicine/sc-rpy:1.0"
    label "process_medium"


    input:
    tuple val(meta), path(adata_umap), val(leiden_resolutions), path(adata_celltypist), path(adata_entropy)

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
    celltypist_adatas = "${adata_celltypist}".split(" ")
    entropy_adatas = "${adata_entropy}".split(" ")

    adata_umap = sc.read_h5ad("${adata_umap}")
    for res, celltypist_path, entropy_path in zip(resolutions, celltypist_adatas, entropy_adatas):
        tmp_adata = sc.read_h5ad(celltypist_path)
        adata_umap.obs[f"leiden_{res:.2f}"] = tmp_adata.obs["leiden"]
        adata_umap.obs[f"leiden_{res:.2f}_celltypist_majority"] = tmp_adata.obs["celltypist_majority"]
        tmp_adata = sc.read_h5ad(entropy_path)
        if "entropy" in tmp_adata.obs.columns:
            adata_umap.obs[f"leiden_{res:.2f}_entropy"] = tmp_adata.obs["entropy"]

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

        ch_resolution = CELLTYPIST_MAJORITY.out.join(
            ENTROPY.out, by: [0, 1]
        ) // meta, resolution, celltypist, entropy

        MERGE_UMAP_LEIDEN(
            UMAP.out.join(
                ch_resolution.groupTuple(by: 0), by: 0
            )
        )

    emit:
        MERGE_UMAP_LEIDEN.out
}
