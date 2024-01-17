include { NEIGHBORS } from "../modules/neighbors.nf"
include { UMAP } from "../modules/umap.nf"
include { LEIDEN } from "../modules/leiden.nf"
include { ENTROPY } from "../modules/entropy.nf"
include { CELLTYPIST_MAJORITY } from "../modules/celltypist.nf"
include { MERGE_CLUSTERING } from "../modules/merge_clustering.nf"


workflow CLUSTERING {
    take:
        ch_adata
        ch_leiden_resolutions
        ch_celltypist

    main:
        NEIGHBORS(ch_adata)
        UMAP(NEIGHBORS.out)
        LEIDEN(NEIGHBORS.out.combine(ch_leiden_resolutions))
        ENTROPY(LEIDEN.out)
        CELLTYPIST_MAJORITY(LEIDEN.out, ch_celltypist)

        ch_resolution = CELLTYPIST_MAJORITY.out.join(
            ENTROPY.out, by: [0, 1]
        ) // meta, resolution, celltypist, entropy

        MERGE_CLUSTERING(
            UMAP.out.join(
                ch_resolution.groupTuple(by: 0), by: 0
            )
        )

    emit:
        MERGE_CLUSTERING.out
}
