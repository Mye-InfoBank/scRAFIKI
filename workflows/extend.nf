// Modules
include { CELLTYPIST } from "../modules/celltypist"
include { CELL_CYCLE } from "../modules/cell_cycle"
include { CELL_QC    } from "../modules/cell_qc"
include { INTEGRATE_SCARCHES } from "../modules/integrate_scarches"

// Workflows
include { PREPROCESSING } from "../subworkflows/preprocessing"
include { COUNTS } from "../subworkflows/counts"
include { DOUBLETS } from "../subworkflows/doublets.nf"

workflow EXTEND {
    if (!params.samplesheet) { 
        exit 1, 'Samplesheet not specified!'
    }

    ch_samplesheet = file(params.samplesheet)

    if (!params.base) {
        exit 1, 'Base path not specified!'
    }

    ch_base = Channel.value(file(params.base)).map{ base -> [[id: "base"], base]}

    if (!params.model) {
        exit 1, 'Model not specified!'
    }

    ch_model = Channel.value(file(params.model)).map{ model -> [[id: "model"], model]}

    PREPROCESSING(ch_samplesheet, ch_base, true)

    ch_adata_intersection = PREPROCESSING.out.intersection
    ch_adata_union = PREPROCESSING.out.union
    ch_adata_transfer = PREPROCESSING.out.transfer
    ch_hvgs = PREPROCESSING.out.hvgs

    COUNTS(ch_adata_union, params.normalization_method)

    CELLTYPIST(
        ch_adata_intersection,
        params.celltypist_model ?: ""
    )

    CELL_CYCLE(
        ch_adata_intersection,
        "human"
    )

    CELL_QC(ch_adata_intersection)

    INTEGRATE_SCARCHES(
        ch_adata_transfer,
        ch_model,
        params.has_celltypes
    )

    DOUBLETS(
        ch_adata_transfer,
        ch_hvgs,
        INTEGRATE_SCARCHES.out.model,
        INTEGRATE_SCARCHES.out.integrated,
        COUNTS.out
    )

    ch_adata = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obs = Channel.empty()
    ch_var = Channel.empty()
    
    emit:
        adata = ch_adata
        obsm = ch_obsm
        obs = ch_obs
        var = ch_var
}