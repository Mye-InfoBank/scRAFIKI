include { CLEAN_ADATA      } from '../modules/clean_adata.nf'
include { SPLIT_CATEGORIES } from '../modules/split_categories.nf'
include { CLUSTERING       } from '../subworkflows/clustering.nf'
include { MERGE            } from '../modules/merge.nf'

workflow SUB {
    ch_adata_input = Channel.fromPath(params.input)
        .map{ adata -> [[id: 'input'], adata]}

    ch_categories = params.category_annotation ? Channel.fromPath(params.category_annotation) : []

    CLEAN_ADATA(ch_adata_input, params.integration)
    ch_adata = CLEAN_ADATA.out.adata

    SPLIT_CATEGORIES(ch_adata, ch_categories, params.category)

    ch_categoric_adata = SPLIT_CATEGORIES.out
        .map{ meta, adatas -> adatas }
        .flatten()
        .map{ adata -> [[id: adata.simpleName], adata]}
    
    CLUSTERING(
        ch_categoric_adata,
        Channel.from(params.leiden_resolutions),
        Channel.empty(),
        Channel.value(params.entropy_initial_smoothness)
    )

    ch_obsm = CLUSTERING.out.obsm
        .mix(CLEAN_ADATA.out.embedding)
        .mix(CLEAN_ADATA.out.umap)
        .map{ meta, obsm -> obsm}
        .collect()

    ch_obs = CLUSTERING.out.obs
        .map{ meta, obs -> obs}
        .collect()

    emit:
        adata = ch_adata
        obsm = ch_obsm
        obs = ch_obs
}