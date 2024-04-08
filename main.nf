#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

// Workflows
include { BUILD  } from "./workflows/build"
include { EXTEND } from "./workflows/extend"
include { SUB    } from "./workflows/sub"

// Modules
include { MERGE               } from "./modules/merge.nf"
include { RENAME_INTEGRATIONS } from "./modules/rename_integrations.nf"

workflow {
    if (params.mode == "build") {
        BUILD()
        subworkflow = BUILD.out
    } else if (params.mode == "extend") {
        EXTEND()
        subworkflow = EXTEND.out
    } else if (params.mode == "sub") {
        SUB()
        subworkflow = SUB.out
    } else {
        error "Invalid mode: ${params.mode}. Please enable one of the following profiles: 'build', 'extend', or 'sub'"
    }


    MERGE (
        subworkflow.adata,
        subworkflow.obsm.map{ meta, obsm -> obsm}.collect(),
        subworkflow.obs.map{ meta, obs -> obs}.collect(),
        subworkflow.var.map{ meta, var -> var}.collect()
    )

    RENAME_INTEGRATIONS(MERGE.out.adata)
}
