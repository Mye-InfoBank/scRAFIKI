// first only scTransform implemented
process NORMALIZE {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy"
    label "process_high"
    label "process_high_memory"

    input:
    tuple val(meta), path(adata) 

    output:
    tuple val(meta), path("normalized.h5ad")

    script:
    """
    normalize.py --input ${adata} --output normalized.h5ad
    """
}
