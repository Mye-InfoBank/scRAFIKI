process IDENTIFY_HVGS {
    tag "$meta.id"

    container = "bigdatainbiomedicine/sc-rapids:0.5"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    val(n_hvgs)
    
    output:
    tuple val(meta), path("${meta.id}.hvgs.h5ad")
    
    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    try:
      import rapids_singlecell as sc
      cuda = True
    except:
      import scanpy as sc
      cuda = False

    adata = ad.read_h5ad("${adata}")

    if cuda:
      sc.utils.anndata_to_GPU(adata)

    sc.pp.highly_variable_genes(adata,
                                n_top_genes=${n_hvgs},
                                flavor="seurat_v3",
                                batch_key="batch")

    adata.write_h5ad("${meta.id}.hvgs.h5ad")
    """
}