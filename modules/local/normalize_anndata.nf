// first only scTransform implemented
process NORMALIZATION {
    container = "bigdatainbiomedicine/sctransform"

    input:
    tuple val(id), val(rep), path(adata) 

    output:
        tuple val(id), val(rep), path("${id}_${rep}_normalized.h5ad")


    script:
    """
    normalization.py ${adata} ${id}_${rep}_normalized.h5ad
    """
}
