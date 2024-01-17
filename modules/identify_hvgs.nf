process IDENTIFY_HVGS {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    val(n_hvgs)
    
    output:
    tuple val(meta), path("${meta.id}.hvgs.h5ad")
    
    script:
    """
    #!/opt/conda/bin/python

    import scanpy as sc

    adata = sc.read_h5ad("${adata}")
    sc.pp.highly_variable_genes(adata,
                                n_top_genes=${n_hvgs},
                                flavor="seurat_v3",
                                batch_key="batch")

    adata.write_h5ad("${meta.id}.hvgs.h5ad")
    """
}