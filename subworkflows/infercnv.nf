include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { H5AD_TO_SCE }  from "../modules/local/scconversion/main.nf"

process RUN_SCEVAN {
    input:
    tuple val(id), path(input_file)

    output:
    path "*", emit: output_file

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    scevan_parallel.R ${input_file} ${task.cpus} ${input_file.baseName} > ${id}.log 2>&1 $ignore_exit_code
    """
}

process RUN_INFERCNVPY {
    input:
    tuple val(id), path(input_file)
    path(gtffile)

    output:
    path "*", emit: output_file

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    infercnvpy_parallel.py ${input_file} ${gtffile} > ${id}.log 2>&1 $ignore_exit_code
    """
}

workflow infercnv {
    take: adata_annotated

    main:
    ch_adata_annotated = adata_annotated.map{ it -> [it.baseName, it]}
    SPLIT_ANNDATA(ch_adata_annotated, "patient")

    ch_adatas_by_patient = SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it]}
    H5AD_TO_SCE(ch_adatas_by_patient)

    RUN_INFERCNVPY(
        ch_adatas_by_patient,
        file("${baseDir}/data/10_references/gencode.v33.primary_assembly.annotation.gtf", checkIfExists: true)
    )
    RUN_SCEVAN(H5AD_TO_SCE.out.sce)
}


