process DEDOUBLET_ADATA {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(solo_annotations)
    
    output:
        tuple val(meta), path("${meta.id}.dedup.h5ad")
    
    script:
    """
    #!/opt/conda/bin/python

    import pandas as pd
    import anndata as ad

    adata = ad.read_h5ad("${adata}")
    solo_paths = "${solo_annotations.join(" ")}".split(" ")

    solo_annotation = pd.concat([pd.read_pickle(path) for path in solo_paths])
    solo_annotation = solo_annotation.reindex(adata.obs_names)

    # Keep only cells with "singlet" in the "doublet_label" column
    adata = adata[solo_annotation["doublet_label"] == "singlet", :]

    # Save the AnnData object
    adata.write_h5ad("${meta.id}.dedup.h5ad")
    """
}