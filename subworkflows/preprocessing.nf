include { check_samplesheet } from '../modules/check_samplesheet'

include { RDS_TO_H5AD } from "../modules/rds_to_h5ad.nf"
include { PREPROCESS } from "../modules/preprocess.nf"
include { COLLECT_PROBLEMS } from "../modules/collect_problems.nf"
include { STOP_IF_PROBLEMS } from "../modules/stop_if_problems.nf"
include { GENES_UPSET } from "../modules/genes_upset.nf"
include { MERGE_DATASETS } from "../modules/merge_datasets.nf"
include { COMPOSITION } from "../modules/composition.nf"
include { DISTRIBUTION } from "../modules/distribution.nf"
include { IDENTIFY_HVGS } from "../modules/identify_hvgs.nf"

workflow PREPROCESSING {
    take:
        ch_samplesheet
        ch_base

    main:
        ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))
            .branch { 
                h5ad: it[2] == "h5ad"
                rds: it[2] == "rds"
            }
        
        RDS_TO_H5AD(ch_samples.rds.map{ meta, rds, format -> [meta, rds]})

        PREPROCESS(
            ch_samples.h5ad.map{ meta, adata, format -> [meta, adata]}.mix(RDS_TO_H5AD.out),
            params.custom_metadata,
        )

        COLLECT_PROBLEMS(PREPROCESS.out.problems.map{ meta, problems -> problems}.collect())
        STOP_IF_PROBLEMS(COLLECT_PROBLEMS.out)

        GENES_UPSET(PREPROCESS.out.adata.mix(ch_base).map{ meta, adata -> adata }.collect())

        MERGE_DATASETS(
            PREPROCESS.out.adata.map{ meta, adata -> adata }.collect(),
            ch_base.collect(),
            params.min_cells,
            params.mode == "build" ? params.custom_hvgs : []
        )
        
        ch_adata_intersection = MERGE_DATASETS.out.intersection
            .map{ adata -> [[id: "intersection"], adata] }

        ch_adata_union = MERGE_DATASETS.out.union
            .map{ adata -> [[id: "union"], adata] }

        ch_adata_transfer = MERGE_DATASETS.out.transfer
            .map{ adata -> [[id: "transfer"], adata] }

        COMPOSITION(ch_adata_intersection)
        DISTRIBUTION(ch_adata_intersection)

        if (params.mode == "build") {
            IDENTIFY_HVGS(
                ch_adata_intersection,
                params.integration_hvgs,
                params.custom_hvgs
            )
            ch_hvgs = IDENTIFY_HVGS.out
        } else {
            ch_hvgs = Channel.empty()
        }



    emit:
        intersection = ch_adata_intersection
        union        = ch_adata_union
        transfer     = ch_adata_transfer
        hvgs         = ch_hvgs
}