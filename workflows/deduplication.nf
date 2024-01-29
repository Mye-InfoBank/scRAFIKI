include { SOLO } from "../modules/solo.nf"
include { DEDUPLICATE_ADATA as DEDUPLICATE_INTEGRATIONS } from "../modules/deduplicate_adata.nf"
include { DEDUPLICATE_ADATA as DEDUPLICATE_RAW } from "../modules/deduplicate_adata.nf"
include { EXTRACT_EMBEDDING } from "../modules/extract_embedding.nf"


workflow DEDUPLICATION {
    take:
        ch_hvgs
        ch_scanvi_model
        ch_integrations
        ch_raw
        ch_batches
    
    main:
        SOLO(
            ch_hvgs,
            ch_scanvi_model,
            ch_batches
        )

        DEDUPLICATE_INTEGRATIONS(
            ch_integrations,
            SOLO.out
        )

        EXTRACT_EMBEDDING(DEDUPLICATE_INTEGRATIONS.out)

        DEDUPLICATE_RAW(
            ch_raw,
            SOLO.out
        )

    emit:
        solo = SOLO.out
        integrations = DEDUPLICATE_INTEGRATIONS.out
        raw = DEDUPLICATE_RAW.out
        obsm = EXTRACT_EMBEDDING.out
}