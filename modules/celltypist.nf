process CELLTYPIST {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_high"
    label "process_high_memory"

    input:
    tuple val(meta), path(adata)
    val(model)
    
    output:
    tuple val(meta), path("${meta.id}.celltypist.pkl")

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    annotate_celltypist.py --input ${adata} --output ${meta.id}.celltypist.pkl --model ${model}
    """
}

process CELLTYPIST_MAJORITY {
    tag "${meta.id}"

    container "bigdatainbiomedicine/sc-rpy:1.0"
    label "process_medium"

    input:
    tuple val(meta), path(adata)
    tuple val(meta2), path(celltypist)
    
    output:
    tuple val(meta), path("${meta.id}.majority.pkl")
    
    script:
    """
    celltypist_majority.py --input_clustering ${adata} --clustering_key ${meta.id} --input_celltypist ${celltypist} --output ${meta.id}.majority.pkl
    """
}