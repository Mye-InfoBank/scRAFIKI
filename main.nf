#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

// Modules
include { CELLTYPIST } from "./modules/celltypist.nf"
include { CELL_CYCLE } from "./modules/cell_cycle.nf"
include { MERGE } from "./modules/merge.nf"

// Workflows
include { PREPROCESSING } from "./workflows/preprocessing.nf"
include { COUNTS } from "./workflows/counts.nf"
include { INTEGRATION } from "./workflows/integration.nf"
include { DEDUPLICATION } from "./workflows/deduplication.nf"
include { CLUSTERING } from "./workflows/clustering.nf"

if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified!' }

workflow {
    PREPROCESSING(ch_samplesheet)

    ch_preprocessed = PREPROCESSING.out.simple
    ch_hvgs = PREPROCESSING.out.hvgs
    ch_batches = PREPROCESSING.out.batches

    COUNTS(ch_preprocessed, params.normalization_method)

    CELLTYPIST(
        ch_preprocessed,
        params.celltypist_model ?: ""
    )

    CELL_CYCLE(
        ch_preprocessed,
        "human"
    )

    INTEGRATION(
        ch_hvgs,
        Channel.from(params.integration_methods).mix(Channel.value("unintegrated")),
        Channel.value(params.benchmark_hvgs)
    )

    DEDUPLICATION(
        ch_hvgs,
        INTEGRATION.out.scanvi_model,
        INTEGRATION.out.integrated,
        ch_preprocessed,
        ch_batches.collect()
    )

    CLUSTERING(
        DEDUPLICATION.out.integrations,
        Channel.from(params.leiden_resolutions),
        CELLTYPIST.out,
        Channel.value(params.entropy_initial_smoothness)
    )

    ch_obs = CLUSTERING.out.obs.mix(
        CELL_CYCLE.out, CELLTYPIST.out, DEDUPLICATION.out.solo
    )

    ch_obsm = CLUSTERING.out.obsm.mix(
        DEDUPLICATION.out.obsm
    )

    MERGE (
        DEDUPLICATION.out.raw,
        COUNTS.out,
        ch_obsm.map{ meta, obsm -> obsm}.collect(),
        ch_obs.map{ meta, obs -> obs}.collect()
    )
}
