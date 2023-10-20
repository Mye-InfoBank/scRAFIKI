// first only scTransform implemented
process NORMALIZATION {

    input:
    tuple val(id), val(rep), path(adata)

    output:
        tuple val(id), val(rep), path("${id}_${rep}_normalized.h5ad")

    script:

    """
    #!/usr/bin/env python

    import scanpy as sc
    import anndata as ad
    from pysctransform import SCTransform

    adata = sc.read_h5ad("${adata}")

    # Get pearson residuals for 3K highly variable genes: default value
    residuals = SCTransform(adata, var_features_n=3000)
    adata.obsm["pearson_residuals"] = residuals

    """

}
