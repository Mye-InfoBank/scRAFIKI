process CELLTYPIST {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy"

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