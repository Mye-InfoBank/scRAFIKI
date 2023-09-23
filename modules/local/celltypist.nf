process CELLTYPIST {
    input:
    tuple val(id), val(rep), path(adata)
    val(perform_majority_voting)
    val(model)
    
    output:
    tuple val(id), val(rep), path("${id}_${rep}_celltypist.h5ad")
    
    script:
    def notebook = "${baseDir}/analyses/30_annotate_scrnaseq_data/31_celltypist.py"
    """
    python3 ${notebook} --input ${adata} --output ${id}_${rep}_celltypist.h5ad --model ${model} --majority_voting ${perform_majority_voting}
    """
}