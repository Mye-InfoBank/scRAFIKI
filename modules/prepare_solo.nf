process PREPARE_SOLO {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    tuple val(meta2), path(adata_core)
    tuple val(meta3), path(hvgs)
    
    output:
    tuple val(meta), path("${meta.id}.hvg.h5ad"), emit: adata
    tuple val(meta), path("batches.txt"), emit: batches
    
    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad("${adata}")
    adata_core = ad.read_h5ad("${adata_core}")
    df_hvgs = pd.read_pickle("${hvgs}")

    known_cell_types = adata_core.obs["cell_type"].unique()
    adata.obs["cell_type"] = adata.obs["cell_type"].map(lambda original: original if original in known_cell_types else "Unknown")

    adata_hvg = adata[:, df_hvgs[df_hvgs["highly_variable"]].index]

    with open("batches.txt", "w") as f:
        f.write("\\n".join(adata_hvg.obs["batch"].unique()))
    
    adata_hvg.write_h5ad("${meta.id}.hvg.h5ad")
    """
}