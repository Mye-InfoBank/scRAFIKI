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
        NEIGHBORS(ch_adata)
        UMAP(NEIGHBORS.out)
        LEIDEN(NEIGHBORS.out.combine(ch_leiden_resolutions))

        ch_leiden_per_integration = LEIDEN.out.table.map{meta, table -> [meta.integration, table]}.groupTuple()
        ch_integations = ch_adata   .map{meta, adata -> [meta.integration, adata]}
                                    .join(ch_leiden_per_integration)
                                    .map{integration, adata, tables -> [[id: integration], adata, tables]}

        SC_HPL(ch_integations)

        ENTROPY(LEIDEN.out.adata, ch_entropy_smoothness)
        CELLTYPIST_MAJORITY(LEIDEN.out.table, ch_celltypist)

        ch_clustering = LEIDEN.out.table.mix(
            ENTROPY.out, CELLTYPIST_MAJORITY.out
        )

    emit:
        obs = ch_clustering
        obsm = UMAP.out
}
