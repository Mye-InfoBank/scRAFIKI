process MERGE_EXTENDED {
    tag "${meta.id}"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
    tuple val(meta), path(extension)
    tuple val(meta2), path(core)
    
    output:
    tuple val(meta), path("${meta.id}.h5ad")

    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    
    adata_core = ad.read_h5ad("${core}")
    adata_extension = ad.read_h5ad("${extension}")

    adata_extended = ad.concat([adata_core, adata_extension])

    adata_extended.write_h5ad("${meta.id}.h5ad")
    """
}