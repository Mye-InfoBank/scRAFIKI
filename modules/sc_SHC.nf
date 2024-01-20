process SCSHC_CLUSTERING {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.1"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), val("${clustering_key}"), path("${meta.id}.${clustering_key}.clustering.h5ad")

    script:
    clustering_key = "scSHC"
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import anndata2ri
    import rpy2
    import rpy2.robjects as ro
    scSHC = ro.packages.importr('scSHC')

    adata = ad.read_h5ad("${adata}")
    sce = anndata2ri.py2rpy(adata)
    expression = ro.r("counts")(sce)
    batches = anndata2ri.py2rpy(adata.obs['batch'])

    n_genes = min(adata.shape[1], 2500)

    result = scSHC.scSHC(expression, batch=batches, cores=${task.cpus}, num_features=n_genes)
    adata.obs["${clustering_key}"] = anndata2ri.rpy2py(result[0])
    adata.obs["${clustering_key}"] = adata.obs["${clustering_key}"].astype(int).astype(str).astype("category")

    adata.write_h5ad("${meta.id}.${clustering_key}.clustering.h5ad")
    """
}

process SCSHC_CLUSTERING_QC {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.1"
    label "process_medium"

    input:
    tuple val(meta), val(clustering_key), path(adata)

    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import anndata2ri
    import rpy2
    import rpy2.robjects as ro
    scSHC = ro.packages.importr('scSHC')

    adata = ad.read_h5ad("${adata}")
    sce = anndata2ri.py2rpy(adata)
    expression = ro.r("counts")(sce)
    batches = anndata2ri.py2rpy(adata.obs['batch'])
    clusters = anndata2ri.py2rpy(adata.obs['${clustering_key}'])

    n_genes = min(adata.shape[1], 2500)

    result = scSHC.testClusters(expression, clusters, batch=batches, cores=${task.cpus}, num_features=n_genes)

    # TODO: What to do with the result?
    """
}