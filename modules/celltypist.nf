process CELLTYPIST {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_high"
    label "process_high_memory"

    input:
    tuple val(meta), path(adata)
    val(model)
    
    output:
    tuple val(meta), path("${meta.id}.celltypist.h5ad")
    
    script:
    """
    annotate_celltypist.py --input ${adata} --output ${meta.id}.celltypist.h5ad --model ${model}
    """
}

process CELLTYPIST_MAJORITY {
    tag "${meta.id}"

    container "bigdatainbiomedicine/sc-rpy"
    label "process_medium"

    input:
    tuple val(meta), val(resolution), path(adata)
    tuple val(meta2), path(celltypist)
    
    output:
    tuple val(meta), val(resolution), path("${meta.id}.res_${resolution}.majority.h5ad")
    
    script:
    """
    celltypist_majority.py --input_clustering ${adata} --input_celltypist ${celltypist} --output ${meta.id}.res_${resolution}.majority.h5ad
    """
}