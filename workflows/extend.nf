// Modules
include { CELLTYPIST } from "../modules/celltypist.nf"
include { CELL_CYCLE } from "../modules/cell_cycle.nf"
include { CELL_QC    } from "../modules/cell_qc.nf"

// Workflows
include { PREPROCESSING } from "../subworkflows/preprocessing.nf"
include { COUNTS } from "../subworkflows/counts.nf"


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

    PREPROCESSING(ch_samplesheet, ch_base)

    ch_adata_intersection = PREPROCESSING.out.intersection
    ch_adata_union = PREPROCESSING.out.union
    ch_adata_transfer = PREPROCESSING.out.transfer

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

    ch_adata = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obs = Channel.empty()
    
    emit:
        adata = ch_adata
        obsm = ch_obsm
        obs = ch_obs
}