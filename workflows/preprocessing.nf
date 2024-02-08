include { check_samplesheet } from '../modules/check_samplesheet'

include { RDS_TO_H5AD } from "../modules/rds_to_h5ad.nf"
include { FILTER } from "../modules/filter.nf"
include { GENES_UPSET } from "../modules/genes_upset.nf"
include { MERGE_DATASETS } from "../modules/merge_datasets.nf"
include { COMPOSITION } from "../modules/composition.nf"
include { DISTRIBUTION } from "../modules/distribution.nf"
include { IDENTIFY_HVGS } from "../modules/identify_hvgs.nf"

workflow PREPROCESSING {
    take:
        ch_samplesheet

    main:
        ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))
            .branch { 
                h5ad: it[2] == "h5ad"
                rds: it[2] == "rds"
            }
        
        RDS_TO_H5AD(ch_samples.rds.map{ meta, rds, format -> [meta, rds]})

        FILTER(ch_samples.h5ad.map{ meta, adata, format -> [meta, adata]}.mix(RDS_TO_H5AD.out))
        GENES_UPSET(FILTER.out.map{ meta, adata -> adata }.collect())
        MERGE_DATASETS(FILTER.out.flatMap{ meta, adata -> adata }.collect())

        ch_preprocessed = MERGE_DATASETS.out.adata
            .map{ adata -> [[id: "preprocessed"], adata] }

        ch_batches = MERGE_DATASETS.out.batches
            .splitText()
            // Remove \n
            .map{ batch -> batch.replace("\n", "") }

        COMPOSITION(ch_preprocessed)
        DISTRIBUTION(ch_preprocessed)

        IDENTIFY_HVGS(
            ch_preprocessed,
            params.integration_hvgs
        )

    emit:
        simple = ch_preprocessed
        hvgs = IDENTIFY_HVGS.out
        batches = ch_batches
}