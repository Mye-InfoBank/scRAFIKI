process GET_BATCHES {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.2"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    
    output:
    tuple val(meta), path("batches.txt")
    
    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad("${adata}")

    with open("batches.txt", "w") as f:
        f.write("\\n".join(adata.obs["batch"].unique()))
    """
}