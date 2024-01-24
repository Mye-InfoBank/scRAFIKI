#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

// Modules
include { CELLTYPIST } from "./modules/celltypist.nf"
include { SOLO } from "./modules/solo.nf"
include { CELL_CYCLE } from "./modules/cell_cycle.nf"
include { MERGE } from "./modules/merge.nf"

// Workflows
include { PREPROCESSING } from "./workflows/preprocessing.nf"
include { COUNTS } from "./workflows/counts.nf"
include { INTEGRATION } from "./workflows/integration.nf"
include { CLUSTERING } from "./workflows/clustering.nf"
include { BENCHMARKING } from "./workflows/benchmarking.nf"

if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified!' }

workflow {
    PREPROCESSING(ch_samplesheet)

    ch_preprocessed = PREPROCESSING.out.simple
    ch_hvgs = PREPROCESSING.out.hvgs

    COUNTS(ch_preprocessed, params.normalization_method)

    CELLTYPIST(
        ch_preprocessed,
        params.celltypist_model ?: ""
    )

    INTEGRATION(
        ch_hvgs,
        Channel.from(params.integration_methods).mix(Channel.value("unintegrated"))
    )

    BENCHMARKING(
        ch_preprocessed,
        INTEGRATION.out.integrated_types,
        params.benchmark_hvgs
    )

    SOLO(
        ch_hvgs,
        INTEGRATION.out.scanvi_model
    )

    CELL_CYCLE(
        ch_preprocessed,
        "human"
    )

    CLUSTERING(
        INTEGRATION.out.integrated,
        Channel.from(params.leiden_resolutions),
        CELLTYPIST.out
    )

    ch_obs = CLUSTERING.out.obs.mix(
        SOLO.out, CELL_CYCLE.out, CELLTYPIST.out
    )

    ch_obsm = CLUSTERING.out.obsm.mix(
        INTEGRATION.out.obsm
    )

    MERGE(
        ch_preprocessed,
        COUNTS.out,
        ch_obsm.map{ meta, obsm -> obsm}.collect(),
        ch_obs.map{ meta, obs -> obs}.collect()
        )

}
