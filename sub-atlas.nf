include { SPLIT_CATEGORIES } from './modules/split_categories.nf'
include { CLUSTERING } from './workflows/clustering.nf'
include { EXTRACT_EMBEDDING } from './modules/extract_embedding.nf'
include { MERGE } from './modules/merge.nf'

ch_adata = Channel.fromPath(params.input)
    .map{ adata -> [[id: 'input'], adata]}
ch_categories = Channel.fromPath(params.category_annotation)

workflow {
    SPLIT_CATEGORIES(ch_adata, ch_categories, params.category)

    ch_categoric_adata = SPLIT_CATEGORIES.out
        .map{ meta, adatas -> adatas }
        .flatten()
        .map{ adata -> [[id: adata.simpleName], adata]}
    
    CLUSTERING(
        ch_categoric_adata,
        Channel.from(params.leiden_resolutions),
        Channel.empty(),
        Channel.value(params.entropy_initial_smoothness),
        Channel.value(params.integration.split("_")[1])
    )

    EXTRACT_EMBEDDING(
        ch_adata,
        "X_" + params.integration
    )

    ch_obsm = CLUSTERING.out.obsm
        .mix(EXTRACT_EMBEDDING.out)
        .map{ meta, obsm -> obsm}
        .collect()

    ch_obs = CLUSTERING.out.obs
        .map{ meta, obs -> obs}
        .collect()

    MERGE(
        ch_adata,
        ch_obsm,
        ch_obs
    )
}