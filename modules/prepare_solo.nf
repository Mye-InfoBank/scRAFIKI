process PREPARE_SOLO {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    tuple val(meta2), path(hvgs)
    
    output:
    tuple val(meta), path("${meta.id}.hvg.h5ad"), emit: adata
    tuple val(meta), path("batches.txt"), emit: batches
    
    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad("${adata}")
    df_hvgs = pd.read_pickle("${hvgs}")

    adata_hvg = adata[:, df_hvgs[df_hvgs["highly_variable"]].index]

    with open("batches.txt", "w") as f:
        f.write("\\n".join(adata_hvg.obs["batch"].unique()))
    
    adata_hvg.write_h5ad("${meta.id}.hvg.h5ad")
    """
}