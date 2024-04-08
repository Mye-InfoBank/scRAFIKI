include { GET_BATCHES } from "../modules/get_batches"
include { SOLO } from "../modules/solo"
include { DEDOUBLET_ADATA as DEDOUBLET_INTEGRATIONS } from "../modules/dedoublet_adata"
include { DEDOUBLET_ADATA as DEDOUBLET_COUNTS } from "../modules/dedoublet_adata"
include { EXTRACT_EMBEDDING } from "../modules/extract_embedding"


workflow DOUBLETS {
    take:
        ch_adata
        ch_model
        ch_integrations
        ch_raw
    
    main:
        GET_BATCHES(
            ch_adata
        )

        ch_batches = GET_BATCHES.out
            .splitText().flatten()
            .map{ batch -> batch.replace("\n", "") }

        SOLO(
            ch_adata,
            ch_model.collect(),
            params.has_celltypes,
            ch_batches
        )

        ch_solo_annotations = SOLO.out
            .map{ meta, annotation -> annotation }.collect()
            .map{ annotations -> [[id: "solo"], annotations]}

        DEDOUBLET_INTEGRATIONS(
            ch_integrations,
            ch_solo_annotations
        )

        EXTRACT_EMBEDDING(DEDOUBLET_INTEGRATIONS.out)

        DEDOUBLET_COUNTS(
            ch_raw,
            ch_solo_annotations
        )

    emit:
        integrations = DEDOUBLET_INTEGRATIONS.out
        counts = DEDOUBLET_COUNTS.out
        obsm = EXTRACT_EMBEDDING.out
}