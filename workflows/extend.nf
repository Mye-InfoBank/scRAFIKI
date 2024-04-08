// Modules
include { CELLTYPIST } from "../modules/celltypist"
include { CELL_CYCLE } from "../modules/cell_cycle"
include { CELL_QC    } from "../modules/cell_qc"
include { INTEGRATE_SCARCHES } from "../modules/integrate_scarches"
include { MERGE_EXTENDED } from "../modules/merge_extended"

// Workflows
include { PREPROCESSING } from "../subworkflows/preprocessing"
include { COUNTS } from "../subworkflows/counts"
include { DOUBLETS } from "../subworkflows/doublets.nf"
include { CLUSTERING } from "../subworkflows/clustering.nf"


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

    INTEGRATE_SCARCHES(
        ch_adata_transfer,
        ch_model,
        params.has_celltypes
    )

    DOUBLETS(
        ch_adata_transfer,
        INTEGRATE_SCARCHES.out.model,
        INTEGRATE_SCARCHES.out.integrated,
        COUNTS.out
    )

    MERGE_EXTENDED(
        DOUBLETS.out.integrations.collect(),
        ch_base,
        params.has_celltypes ? "scanvi" : "scvi"
    )

    CLUSTERING(
        MERGE_EXTENDED.out,
        Channel.from(params.leiden_resolutions),
        CELLTYPIST.out,
        Channel.value(params.entropy_initial_smoothness)
    )

    ch_obs = CLUSTERING.out.obs.mix(
        CELL_CYCLE.out,
        CELLTYPIST.out,
        CELL_QC.out
    )

    ch_obsm = CLUSTERING.out.obsm.mix(
        DOUBLETS.out.obsm
    )
    
    emit:
        adata = DOUBLETS.out.counts
        obsm = ch_obsm
        obs = ch_obs
}