// first only scTransform implemented
process NORMALIZE {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_high"
    label "process_high_memory"

    input:
    tuple val(meta), path(adata)
    val(normalization)

    output:
    tuple val(meta), path("normalized.h5ad")

    script:
    if (normalization == "scTransform")
        """
        sctransform.py --input ${adata} --output normalized.h5ad
        """
    else if (normalization == "scanpy")
        """
        #!/opt/conda/bin/python

        import scanpy as sc

        adata = sc.read_h5ad("${adata}")

        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)

        adata.write_h5ad("normalized.h5ad")
        """
    else
        """
        echo "Unknown normalization method: ${normalization}"
        echo "Valid options are: scTransform, scanpy"
        exit 1
        """
}
