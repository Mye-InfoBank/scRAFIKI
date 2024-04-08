include { SC_HPL } from "../modules/scHPL"
include { NEIGHBORS } from "../modules/neighbors"
include { UMAP } from "../modules/umap"
include { LEIDEN } from "../modules/leiden"
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

        ENTROPY(LEIDEN.out.adata, ch_entropy_smoothness)
        CELLTYPIST_MAJORITY(LEIDEN.out.table, ch_celltypist)

        ch_clustering = LEIDEN.out.table.mix(
            ENTROPY.out, CELLTYPIST_MAJORITY.out
        )

    emit:
        obs = ch_clustering
        obsm = UMAP.out
}
