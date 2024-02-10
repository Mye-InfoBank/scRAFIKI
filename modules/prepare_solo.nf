process PREPARE_SOLO {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    tuple val(meta2), path(hvgs)
    
    output:
    tuple val(meta), path("*.batch.h5ad")
    
    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad("${adata}")
    df_hvgs = pd.read_pickle("${hvgs}")

    adata_hvg = adata[:, df_hvgs[df_hvgs["highly_variable"]].index]

    for batch in adata_hvg.obs["batch"].unique():
        adata_batch = adata_hvg[adata_hvg.obs["batch"] == batch]
        adata_batch.write_h5ad(f"{batch}.batch.h5ad")
    """
}