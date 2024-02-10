#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

// Modules
include { CELLTYPIST } from "./modules/celltypist.nf"
include { CELL_CYCLE } from "./modules/cell_cycle.nf"
include { CELL_QC    } from "./modules/cell_qc.nf"
include { MERGE } from "./modules/merge.nf"

// Workflows
include { PREPROCESSING } from "./workflows/preprocessing.nf"
include { COUNTS } from "./workflows/counts.nf"
include { INTEGRATION } from "./workflows/integration.nf"
include { DOUBLETS } from "./workflows/doublets.nf"
include { CLUSTERING } from "./workflows/clustering.nf"

if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified!' }

workflow {
    PREPROCESSING(ch_samplesheet)

    ch_adata_integration = PREPROCESSING.out.integration
    ch_adata_intersection = PREPROCESSING.out.intersection
    ch_adata_counts = PREPROCESSING.out.counts
    ch_hvgs = PREPROCESSING.out.hvgs
    ch_batches = PREPROCESSING.out.batches

    COUNTS(ch_adata_counts, params.normalization_method)

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
        ch_hvgs,
        Channel.from(params.integration_methods).mix(Channel.value("unintegrated")),
        Channel.value(params.benchmark_hvgs)
    )

    DOUBLETS(
        ch_hvgs,
        INTEGRATION.out.scanvi_model,
        INTEGRATION.out.integrated,
        COUNTS.out,
        ch_batches.collect()
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
        DOUBLETS.out.solo,
        INTEGRATION.out.scanvi_labels,
        CELL_QC.out
    )

    ch_obsm = CLUSTERING.out.obsm.mix(
        DOUBLETS.out.obsm
    )

    MERGE (
        DOUBLETS.out.counts,
        ch_obsm.map{ meta, obsm -> obsm}.collect(),
        ch_obs.map{ meta, obs -> obs}.collect()
    )
}
