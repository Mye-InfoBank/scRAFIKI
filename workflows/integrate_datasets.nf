include { check_samplesheet } from '../modules/local/check_samplesheet'

include { FILTER } from "../modules/local/filter.nf"


include { SOLO } from "../modules/local/solo/main.nf"
include { NEIGHBORS_LEIDEN_UMAP } from "./neighbors_leiden_umap.nf"
include { MERGE_INTEGRATIONS } from "../modules/local/merge_integrations.nf"
include { BENCHMARK_INTEGRATIONS } from "../modules/local/scIB.nf"
include { DECONTX } from "../modules/local/decontX.nf"
include { CONCAT_ADATA as CONCAT_BATCHES } from "../modules/local/concat_anndata.nf"
include { FILTER_ANNDATA as SPLIT_BATCHES } from "../modules/local/scconversion/main.nf"
include { FILTER_SOLO } from "../modules/local/solo/main.nf"
include { ADATA_METRICS as AFTER_QC_ADATA_METRICS } from "../modules/local/adata_metrics.nf"
include { ADATA_METRICS as FILTERED_ADATA_METRICS } from "../modules/local/adata_metrics.nf"
include { ADATA_METRICS as RAW_ADATA_METRICS } from "../modules/local/adata_metrics.nf"
include { COMBINE_ADATA_METRICS } from "../modules/local/adata_metrics.nf"

include { INTEGRATE } from "../modules/local/integrate.nf"
include { MERGE_DATASETS } from "../modules/local/merge_datasets.nf"

if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified!' }


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

    /*

    ch_adata_merged = MERGE_ALL.out.artifacts.flatten().filter{
        it.getExtension() == "h5ad"
    }

    AFTER_QC_ADATA_METRICS(ch_adata_merged.map{ [[id: "after_qc"], it] })

    ch_adata_merged_meta = ch_adata_merged.collect().map{
        out -> ["all", out]
    }
    */

    ch_integration_methods = Channel.from(params.integration_methods)

    INTEGRATE(
        MERGE_DATASETS.out,
        ch_integration_methods
    )

    /*
  
    ch_scvi_hvg = SCVI.out.adata
    ch_scvi_hvg_model = SCVI.out.scvi_model

    ch_scanvi_hvg = SCANVI.out.adata
    ch_scanvi_hvg_model = SCANVI.out.scvi_model

    ch_batches = MERGE_ALL.out.artifacts.flatten().filter{
            it -> it.getName() == "obs_all.csv"
        }.splitCsv(header : true).map{ it -> it["batch"] ?: "batch" }

    ch_scvi_hvg_complete = ch_scvi_hvg.map{ it[1] }.merge(ch_scvi_hvg_model.map{ it[1] })
        .map{ adata, model -> ["scvi_hvg", "scVI", adata, model, "embed"] }

    ch_scanvi_hvg_complete = ch_scanvi_hvg.map{ it[1] }.merge(ch_scanvi_hvg_model.map{ it[1] })
        .map{ adata, model -> ["scanvi_hvg", "scANVI", adata, model, "embed"] }

    ch_harmony_complete = ch_harmony_adata.map{ adata -> ["harmony", "pca_harmony", adata, null, "embed"] }
    ch_mnn_complete = ch_mnn.map{ adata -> ["mnn", "MNN", adata, null, "embed"] }


    ch_integrations = Channel.empty().mix(
        ch_scvi_hvg_complete,
        ch_scanvi_hvg_complete,
        ch_harmony_complete,
        ch_mnn_complete
    )



    BENCHMARK_INTEGRATIONS(
        ch_adata_merged.combine(ch_integrations.map { [it[0], "X_" + it[1], it[2], it[4]] }),
        params.organism,
        params.scib_fast
    )

    MERGE_INTEGRATIONS(
        ch_adata_merged,
        ch_integrations.map { it[2] }.collect(),
        ch_integrations.map { it[1] }.collect()
    )

    SPLIT_BATCHES(
        ch_batches.combine(MERGE_INTEGRATIONS.out),
        "lambda x: x['batch'] == process_id"
    )

    SOLO(
        ch_scanvi_hvg,
        ch_scanvi_hvg_model,
        ch_batches
    )


    DECONTX(
        SPLIT_BATCHES.out.adata
    )

    FILTER_SOLO(
        DECONTX.out.join(SOLO.out.map{ [it[0], it[1]] }),
    )

    CONCAT_BATCHES(
        FILTER_SOLO.out.map{ it[1] }.collect()
    )

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
