include { check_samplesheet } from '../modules/local/check_samplesheet'

include { FILTER } from "../modules/local/filter.nf"
include { MERGE_DATASETS } from "../modules/local/merge_datasets.nf"

workflow PREPROCESSING {
    take:
        ch_samplesheet

    main:
        ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))
        FILTER(ch_samples)
        MERGE_DATASETS(FILTER.out.flatMap{ meta, adata -> adata }.collect())

        ch_preprocessed = MERGE_DATASETS.out
            .map{ adata -> [[id: "preprocessed"], adata] }

    emit:
        ch_preprocessed
}