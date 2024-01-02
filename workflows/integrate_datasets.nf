include { check_samplesheet } from '../modules/local/check_samplesheet'

include { FILTER } from "../modules/local/filter.nf"


include { SOLO } from "../modules/local/solo/main.nf"
include { NEIGHBORS_LEIDEN_UMAP } from "./neighbors_leiden_umap.nf"
include { MERGE_INTEGRATIONS } from "../modules/local/merge_integrations.nf"
include { BENCHMARK_INTEGRATIONS } from "../modules/local/scIB.nf"
include { DECONTX } from "../modules/local/decontX.nf"
include { CONCAT_ADATA as CONCAT_DECONTX } from "../modules/local/concat_anndata.nf"
include { FILTER_SOLO } from "../modules/local/solo/main.nf"
include { ADATA_METRICS as AFTER_QC_ADATA_METRICS } from "../modules/local/adata_metrics.nf"
include { ADATA_METRICS as FILTERED_ADATA_METRICS } from "../modules/local/adata_metrics.nf"
include { ADATA_METRICS as RAW_ADATA_METRICS } from "../modules/local/adata_metrics.nf"
include { COMBINE_ADATA_METRICS } from "../modules/local/adata_metrics.nf"

include { INTEGRATE } from "../modules/local/integrate.nf"
include { MERGE_DATASETS } from "../modules/local/merge_datasets.nf"
include { SPLIT_BATCHES } from "../modules/local/split_batches.nf"
include { INTEGRATE as INTEGRATE_SCVI } from "../modules/local/integrate.nf"
include { INTEGRATE_SCANVI } from "../modules/local/integrate_scanvi.nf"
include { DROPLET_CORRECT } from "../modules/local/droplet_correct.nf"

if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified!' }

integration_types = [
    "bbknn": "knn",
    "combat": "full",
    "desc": "embed",
    "harmony": "embed",
    "mnn": "full",
    "scanorama": "embed",
    "scanvi": "embed",
    "scvi": "embed",
    "trvaep": "embed",
    "scgen": "full",
]

/**
 * Integrate individual datasets into a single-cell atlas
 *   - Perfom QC on individual datasets
 *   - Perform manual "seed" annotation of two datasets (one SS2, one 10x)
 *   - Perform data integration using SCANVI
 *   - Call doublets using solo
 */
workflow integrate_datasets {


    main:

    ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))
    RAW_ADATA_METRICS(ch_samples)

    FILTER(ch_samples)

    FILTERED_ADATA_METRICS(FILTER.out)

    // MERGE and INTEGRATE all datasets
    MERGE_DATASETS(FILTER.out.flatMap{ meta, adata -> adata }.collect())

    ch_unintegrated = MERGE_DATASETS.out
        .map{ adata -> [[id: "unintegrated"], adata] }

    SPLIT_BATCHES(ch_unintegrated)
    ch_unintegrated_batches = SPLIT_BATCHES.out
        .transpose()
        .map{ meta, adata -> 
            batch = adata.simpleName;
            return [[id: batch], adata]
        }
    
    DECONTX(
        ch_unintegrated_batches
    )

    CONCAT_DECONTX(
        DECONTX.out.map{ it[1] }.collect().map{ [[id: "decontX"], it] }
    )

    ch_integration_methods = Channel.from(params.integration_methods)
        .filter{ it != "scvi" && it != "scanvi" }

    INTEGRATE(
        ch_unintegrated,
        ch_integration_methods
    )

    INTEGRATE_SCVI(
        ch_unintegrated,
        "scvi"
    )

    INTEGRATE_SCANVI(
        ch_unintegrated,
        INTEGRATE_SCVI.out.model
    )

    ch_integrated = INTEGRATE.out.integrated
        .mix(INTEGRATE_SCVI.out.integrated)
        .mix(INTEGRATE_SCANVI.out.integrated)
        .map{ meta, adata -> [meta, adata, integration_types[meta.integration]] }


    if (params.benchmark) {
        BENCHMARK_INTEGRATIONS(
            ch_unintegrated,
            ch_integrated,
            params.organism,
            params.benchmark_hvgs
        )
    }

    MERGE_INTEGRATIONS(
        ch_unintegrated,
        ch_integrated.map { meta, adata, type -> meta.integration }.collect(),
        ch_integrated.map { meta, adata, type -> type }.collect(),
        ch_integrated.map { meta, adata, type -> adata }.collect(),
    )

    ch_merged_integrations = MERGE_INTEGRATIONS.out
        .map{ adata -> [[id: "integrated"], adata] }

    SOLO(
        ch_unintegrated,
        INTEGRATE_SCANVI.out.model
    )

    DROPLET_CORRECT (
        ch_merged_integrations,
        SOLO.out,
        CONCAT_DECONTX.out
    )

    /*
    FILTERED_ADATA_METRICS(
        CONCAT_BATCHES.out.map{ [[id: "no_doublets"], it] }
    )

    COMBINE_ADATA_METRICS(
        AFTER_QC_ADATA_METRICS.out.mix(FILTERED_ADATA_METRICS.out, RAW_ADATA_METRICS.out).collect()
    )

    NEIGHBORS_LEIDEN_UMAP(
        CONCAT_BATCHES.out.map{ ["all", it]},
        ch_integrations.map{ "X_" + it[1] },
        Channel.from(params.clustering_resolutions)
    )

    emit:
        adata_integrated = NEIGHBORS_LEIDEN_UMAP.out.adata

    */

}
