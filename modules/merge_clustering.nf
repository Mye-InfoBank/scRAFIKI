process MERGE_CLUSTERING {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy"
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