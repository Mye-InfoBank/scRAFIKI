process DEDUPLICATE_ADATA {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(solo)
    
    output:
        tuple val(meta), path("${meta.id}.dedup.h5ad")
    
    script:
    """
    #!/opt/conda/bin/python

    import pandas as pd
    import anndata as ad

    adata = ad.read_h5ad("${adata}")
    solo = pd.read_pickle("${solo}")

    # Keep only cells with "singlet" in the "doublet_label" column
    adata = adata[solo["doublet_label"] == "singlet", :]

    # Save the AnnData object
    adata.write_h5ad("${meta.id}.dedup.h5ad")
    """
}