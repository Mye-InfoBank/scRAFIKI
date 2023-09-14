#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

include { integrate_datasets } from "./workflows/integrate_datasets.nf"
include { annotate_dataset } from "./workflows/annotate_dataset.nf"
include { add_additional_datasets } from "./workflows/add_additional_datasets.nf"

workflow {

    integrate_datasets()

    annotate_dataset(integrate_datasets.out.adata_integrated)

    /*
    add_additional_datasets(
        annotate_dataset.out.final_atlas,
        annotate_dataset.out.scanvi_h5ad,
        annotate_dataset.out.scanvi_model
    )
    */
}
