include { check_samplesheet } from '../modules/local/check_samplesheet'

include { SCQC } from "../modules/local/scqc/main"
include { SCQC_MERGE_STATS } from "../modules/local/scqc_merge_stats/main.nf"

include { JUPYTERNOTEBOOK as MERGE_ALL } from "../modules/local/jupyternotebook/main.nf"
include { SCVI } from "../modules/local/scvi/main.nf"
include { SCANVI } from "../modules/local/scvi/main.nf"
include { SOLO } from "../modules/local/solo/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_DOUBLET } from "./neighbors_leiden_umap.nf"
include { JUPYTERNOTEBOOK as MERGE_SOLO }  from "../modules/local/jupyternotebook/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_NODOUBLET } from "./neighbors_leiden_umap.nf"
include { JUPYTERNOTEBOOK as HARMONY }  from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as MNN }  from "../modules/local/jupyternotebook/main.nf"
include { MERGE_INTEGRATIONS } from "../modules/local/merge_integrations.nf"
include { BENCHMARK_INTEGRATIONS } from "../modules/local/scIB.nf"
include { DECONTX } from "../modules/local/decontX.nf"
include { CONCAT_ADATA as CONCAT_BATCHES } from "../modules/local/concat_anndata.nf"
include { FILTER_ANNDATA as SPLIT_BATCHES } from "../modules/local/scconversion/main.nf"
include { FILTER_SOLO } from "../modules/local/solo/main.nf"


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

    SCQC(
        [
            file("${baseDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${baseDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())


    // MERGE and INTEGRATE all datasets
    MERGE_ALL(
        Channel.value([
            [id: "21_merge_all"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/21_merge_all.py")
        ]),
        [
            samplesheet: ch_samplesheet.toString(),
            gene_symbol_table: "gene_symbol_dict.csv"
        ],
        SCQC.out.adata.flatMap{ id, adata -> adata }.mix(
            Channel.value(ch_samplesheet),
            Channel.fromPath("${baseDir}/tables/gene_symbol_dict.csv")
        ).collect()
    )

    ch_adata_merged = MERGE_ALL.out.artifacts.flatten().filter{
        it.getExtension() == "h5ad"
    }


    ch_adata_merged_meta = ch_adata_merged.collect().map{
        out -> ["all", out]
    }

    SCVI(
        ch_adata_merged_meta,
        1, // 1 = use HVG
        ["batch", "dataset", null]
    )

    SCANVI(
        SCVI.out.adata.join(SCVI.out.scvi_model),
        "batch",
        "cell_type"
    )

    HARMONY(
            Channel.value([
            [id: "22_harmony"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/22_harmony.py")
        ]),
        [
        ],
        ch_adata_merged
    )

    ch_harmony_adata = HARMONY.out.artifacts.filter{
        it.getExtension() == "h5ad"
    }

    // ch_harmony_adata.view()

    // Takes very long according to this https://github.com/chriscainx/mnnpy#speed
    if (params.mnn) {
        MNN(
            Channel.value([
                [id: "23_mnn"],
                file("${baseDir}/analyses/20_integrate_scrnaseq_data/23_mnn.py")
            ]),
            [
            ],
            ch_adata_merged
        )

        ch_mnn_adata = MNN.out.artifacts.filter{
            it.getExtension() == "h5ad"
        }

        ch_mnn = Channel.from([
            ["mnn", "MNN", ch_mnn_adata, null]
        ])
    } else {
        ch_mnn = Channel.empty()
    }
  
    ch_scvi_hvg = SCVI.out.adata
    ch_scvi_hvg_model = SCVI.out.scvi_model

    ch_scanvi_hvg = SCANVI.out.adata
    ch_scanvi_hvg_model = SCANVI.out.scvi_model

    ch_batches = MERGE_ALL.out.artifacts.flatten().filter{
            it -> it.getName() == "obs_all.csv"
        }.splitCsv(header : true).filter{
            it -> it["run_solo"] == "True"
        }.map{ it -> it["batch"] }

    ch_scvi_hvg_complete = ch_scvi_hvg.map{ it[1] }.merge(ch_scvi_hvg_model.map{ it[1] })
        .map{ adata, model -> ["scvi_hvg", "scVI", adata, model] }

    ch_scanvi_hvg_complete = ch_scanvi_hvg.map{ it[1] }.merge(ch_scanvi_hvg_model.map{ it[1] })
        .map{ adata, model -> ["scanvi_hvg", "scANVI", adata, model] }

    ch_harmony_complete = ch_harmony_adata.map{ adata -> ["harmony", "pca_harmony", adata, null] }
    ch_mnn_complete = ch_mnn.map{ adata -> ["mnn", "MNN", adata, null] }


    ch_integrations = Channel.empty().mix(
        ch_scvi_hvg_complete,
        ch_scanvi_hvg_complete,
        ch_harmony_complete,
        ch_mnn_complete
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


    /*
    BENCHMARK_INTEGRATIONS(
        MERGE_INTEGRATIONS.out,
        ch_integrations.map { it[1] }.collect()
    )
    */


    DECONTX(
        SPLIT_BATCHES.out.adata
    )

    FILTER_SOLO(
        DECONTX.out.join(SOLO.out.map{ [it[0], it[1]] }),
    )

    CONCAT_BATCHES(
        FILTER_SOLO.out.map{ it[1] }.collect()
    )


    NEIGHBORS_LEIDEN_UMAP_NODOUBLET(
        CONCAT_BATCHES.out,
        "X_scANVI",
        1.0
    )

    emit:
        adata_integrated = NEIGHBORS_LEIDEN_UMAP_NODOUBLET.out.adata.map{ meta, ad -> ad }

}
