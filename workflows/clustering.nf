include { NEIGHBORS } from "../modules/neighbors.nf"
include { UMAP } from "../modules/umap.nf"
include { LEIDEN } from "../modules/leiden.nf"
include { SCSHC_CLUSTERING } from "../modules/sc_SHC.nf"
include { SCSHC_CLUSTERING_QC } from "../modules/sc_SHC.nf"
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
        SCSHC_CLUSTERING(ch_adata.filter{ meta, adata -> meta.integration == "unintegrated" })

        ch_clusterings = LEIDEN.out.mix(SCSHC_CLUSTERING.out)

        SCSHC_CLUSTERING_QC(ch_clusterings)

        ENTROPY(ch_clusterings)
        CELLTYPIST_MAJORITY(ch_clusterings, ch_celltypist)

        ch_cluster_results = CELLTYPIST_MAJORITY.out.join(
            ENTROPY.out, by: [0, 1]
        ) // meta, clustering_key, celltypist, entropy

        ch_clustering_keys = ch_cluster_results.map{ it[1] }.unique()

        MERGE_CLUSTERING(
            UMAP.out.combine(
                ch_cluster_results.groupTuple(by: 0), by: 0
            )
        )

    emit:
        results = MERGE_CLUSTERING.out
        keys = ch_clustering_keys
}
