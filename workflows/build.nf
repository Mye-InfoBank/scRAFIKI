// Modules
include { CELLTYPIST } from "../modules/celltypist.nf"
include { CELL_CYCLE } from "../modules/cell_cycle.nf"
include { CELL_QC    } from "../modules/cell_qc.nf"

// Workflows
include { PREPROCESSING } from "../subworkflows/preprocessing.nf"
include { COUNTS } from "../subworkflows/counts.nf"
include { INTEGRATION } from "../subworkflows/integration.nf"
include { DOUBLETS } from "../subworkflows/doublets.nf"
include { CLUSTERING } from "../subworkflows/clustering.nf"

workflow BUILD {
    if (!params.samplesheet) { 
        exit 1, 'Samplesheet not specified!' 
    }

    ch_samplesheet = file(params.samplesheet) 

    PREPROCESSING(ch_samplesheet, Channel.value([[], []]))

    ch_adata_intersection = PREPROCESSING.out.intersection
    ch_adata_union = PREPROCESSING.out.union
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

    INTEGRATION(
        ch_adata_intersection,
        ch_hvgs,
        Channel.from(params.integration_methods).mix(Channel.value("unintegrated")),
        Channel.value(params.benchmark_hvgs)
    )

    DOUBLETS(
        ch_adata_intersection,
        INTEGRATION.out.model,
        INTEGRATION.out.integrated,
        COUNTS.out
    )

    CLUSTERING(
        DOUBLETS.out.integrations,
        Channel.from(params.leiden_resolutions),
        CELLTYPIST.out,
        Channel.value(params.entropy_initial_smoothness)
    )

    ch_obs = CLUSTERING.out.obs.mix(
        CELL_CYCLE.out,
        CELLTYPIST.out,
        INTEGRATION.out.scanvi_labels,
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
