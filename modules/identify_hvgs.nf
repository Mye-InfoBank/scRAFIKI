process IDENTIFY_HVGS {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    val(n_hvgs)
    val(custom_hvgs)
    
    output:
    tuple val(meta), path("${meta.id}.hvgs.pkl")
    
    script:
    """
    #!/opt/conda/bin/python

    import scanpy as sc

    adata = sc.read_h5ad("${adata}")

    span = 0.3 # default
    worked = False

    while not worked and span <= 1:
        try:
            sc.pp.highly_variable_genes(adata,
                                        n_top_genes=10000,
                                        flavor="seurat_v3",
                                        span=span,
                                        batch_key="batch")
            worked = True
        except:
            span += 0.1
            print(f"Increased span to {span}")
    
    custom_hvgs = "${custom_hvgs.join(" ")}".split(" ")
    custom_hvgs = [x.upper().replace("_", "-").replace(".", "-") for x in custom_hvgs]
    custom_hvgs = [x for x in custom_hvgs if x in adata.var_names]

    if len(custom_hvgs) > 0:
        # Set "highly_variable" to True for custom HVGs
        adata.var.loc[custom_hvgs, "highly_variable"] = True

    adata.var[["highly_variable"]].to_pickle("${meta.id}.hvgs.pkl")
    """
}