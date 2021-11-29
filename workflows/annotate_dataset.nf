
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_COARSE }  from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_FINE }  from "../modules/local/jupyternotebook/main.nf"
include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_CELL_TYPES } from "../subworkflows/neighbors_leiden_umap/main.nf"

// include { PREPARE_ANNDATA as PREPARE_ANNDATA_DE_EPI }  from "../modules/local/scde/main.nf"
// include { MAKE_PSEUDOBULK as MAKE_PSEUDOBULK_EPI }  from "../modules/local/scde/main.nf"
// include { DE_EDGER as DE_EDGER_EPI } from "../modules/local/scde/main.nf"
// include { DE_EDGER as DE_EDGER_EPI_N_CELLS } from "../modules/local/scde/main.nf"
// include { H5AD_TO_SCE as H5AD_TO_SCE_EPI } from "../modules/local/scconversion/main.nf"
// include { DE_MAST_MIXED_EFFECTS as DE_MAST_MIXED_EFFECTS_EPI } from "../modules/local/scde/main.nf"
// include { DE_DREAM as DE_DREAM_EPI } from "../modules/local/scde/main.nf"
// include { DE_DREAM as DE_DREAM_EPI_N_CELLS } from "../modules/local/scde/main.nf"

/**
 * Annotate cell-types of the lung cancer atlas.
 *   - perform a coarse-grained annotation
 *   - perform DE analysis on some clusters
 *   - manually annotate sub-clusters to obtain a fine-grained cell-type annotation.
 */
workflow annotate_dataset {
    take:
        adata_integrated

    main:
    ANNOTATE_CELL_TYPES_COARSE(
        Channel.value([
            [id: "annotate_cell_types_coarse"],
            file("${baseDir}/analyses/30_annotate_scrnaseq_data/31_annotate_cell_types_coarse.py")
        ]),
        [:],
        adata_integrated
    )
    ch_adata_annotated = ANNOTATE_CELL_TYPES_COARSE.out.artifacts
    SPLIT_ANNDATA(
        ch_adata_annotated.map{ it -> [it.baseName, it]},
        "cell_type"
    )
    NEIGHBORS_LEIDEN_UMAP_CELL_TYPES(
        SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it] },
        "X_scANVI",
        Channel.from(0.5, 0.75, 1.0, 1.5)
    )
    ANNOTATE_CELL_TYPES_FINE(
        Channel.value([
            [id: "annotate_cell_types_fine"],
            file("${baseDir}/analyses/30_annotate_scrnaseq_data/32_annotate_cell_types_fine.py")
        ]),
        [
            "input_dir": '.',
            "main_adata": 'adata_cell_type_coarse.h5ad'
        ],
        NEIGHBORS_LEIDEN_UMAP_CELL_TYPES.out.adata.map{ id, adata -> adata }.mix(
            ch_adata_annotated
        ).collect()
    )

    // PREPARE_CELLXGENE(
    //     Channel.value([
    //         [id: "zz_prepare_cellxgene"],
    //         file("${baseDir}/analyses/zz_cellxgene/stats_and_cellxgene.py")
    //     ]),
    //     ["adata_in": "adata_annotated_fine.h5ad"],
    //     ANNOTATE_CELL_TYPES_FINE.out.artifacts
    // )

    // ch_adata_annotated_fine = ANNOTATE_CELL_TYPES_FINE.out.artifacts.flatten().filter(
    //     it -> it.name.contains("h5ad")
    // ).map{ it -> [it.baseName, it]}

    // SPLIT_ANNDATA(ch_adata_annotated_fine, "dataset")
    // H5AD_TO_SEURAT(ch_adata_annotated_fine)

    // ch_epithelial = NEIGHBORS_LEIDEN_UMAP_CELL_TYPES.out.adata.filter{
    //     id, adata -> id.contains("epithelial")
    // }
    // PREPARE_ANNDATA_DE_EPI(
    //     ch_epithelial,
    //     "X",
    //     "leiden_0.50",
    //     ["all", "rest"]
    // )
    // ch_adata_for_de = PREPARE_ANNDATA_DE_EPI.out.adata.flatMap{ id, adatas -> adatas }.map{ it -> [it.baseName, it]}

    // MAKE_PSEUDOBULK_EPI(
    //     ch_adata_for_de,
    //     "patient",
    //     "leiden_0.50",
    //     [10, false]
    // )
    // ch_pseudobulk_leiden = MAKE_PSEUDOBULK_EPI.out.pseudobulk.filter{
    //     id, counts, samplesheet -> samplesheet.text.count("\n") >= 6
    // }
    // // DE_DREAM_EPI(
    // //     ch_pseudobulk_leiden,
    // //     "leiden_0.50",
    // //     "+ (1 | dataset) + (1 | patient:dataset)"
    // // )
    // // DE_DREAM_EPI_N_CELLS(
    // //     ch_pseudobulk_leiden,
    // //     "leiden_0.50",
    // //     "+ n_cells + (1 | dataset) + (1 | patient:dataset)"
    // // )
    // DE_EDGER_EPI(
    //     ch_pseudobulk_leiden,
    //     "leiden_0.50",
    //     " + patient"
    // )
    // DE_EDGER_EPI_N_CELLS(
    //     ch_pseudobulk_leiden,
    //     "leiden_0.50",
    //     " + patient + n_cells"
    // )
    // // H5AD_TO_SCE_EPI(ch_adata_for_de)
    // // DE_MAST_MIXED_EFFECTS_EPI(
    // //     H5AD_TO_SCE_EPI.out.sce,
    // //     "leiden_0.50",
    // //     " + n_genes_by_counts + (1 | patient:dataset) + (1 | dataset)"
    // // )


    // // emit:
    // //     adata_annotated_by_cell_type = ch_adata_annotated_by_cell_type
    // //     adata_annotated_cell_type_coarse = ch_adata_annotated
}
