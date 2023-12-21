#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

include { integrate_datasets } from "./workflows/integrate_datasets.nf"
include { finalize_dataset } from "./workflows/finalize_dataset.nf"
// include { add_additional_datasets } from "./workflows/add_additional_datasets.nf"

workflow {

    integrate_datasets()

    finalize_dataset(integrate_datasets.out.adata_integrated)

    /*
    add_additional_datasets(
        finalize_dataset.out.final_atlas,
        finalize_dataset.out.scanvi_h5ad,
        finalize_dataset.out.scanvi_model
    )
    */
}
