include { SC_HPL } from "../modules/scHPL"
include { NEIGHBORS } from "../modules/neighbors"
include { UMAP } from "../modules/umap"
include { LEIDEN } from "../modules/leiden"
include { SCSHC_CLUSTERING } from "../modules/sc_SHC"
include { SCSHC_CLUSTERING_QC } from "../modules/sc_SHC"
include { ENTROPY } from "../modules/entropy"
include { CELLTYPIST_MAJORITY } from "../modules/celltypist"
include { MERGE_CLUSTERING } from "../modules/merge_clustering"


workflow CLUSTERING {
    take:
        ch_adata
        ch_leiden_resolutions
        ch_celltypist
        ch_entropy_smoothness

    main:
        SC_HPL(ch_adata, [[], []], [[], []])
        NEIGHBORS(ch_adata)
        UMAP(NEIGHBORS.out)
        LEIDEN(NEIGHBORS.out.combine(ch_leiden_resolutions))
        SCSHC_CLUSTERING(ch_adata.filter{ meta, adata -> meta.integration == "unintegrated" })

        ch_clustering_adatas = LEIDEN.out.adata.mix(SCSHC_CLUSTERING.out.adata)
        ch_clustering_tables = LEIDEN.out.table.mix(SCSHC_CLUSTERING.out.table)

        SCSHC_CLUSTERING_QC(ch_clustering_adatas)
        ENTROPY(ch_clustering_adatas, ch_entropy_smoothness)
        CELLTYPIST_MAJORITY(ch_clustering_tables, ch_celltypist)

        ch_clustering = ch_clustering_tables.mix(
            SCSHC_CLUSTERING_QC.out.table, ENTROPY.out, CELLTYPIST_MAJORITY.out
        )

    emit:
        obs = ch_clustering
        obsm = UMAP.out
}
