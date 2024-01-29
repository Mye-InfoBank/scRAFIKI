include { SOLO } from "../modules/solo.nf"
include { DEDOUBLET_ADATA as DEDOUBLET_INTEGRATIONS } from "../modules/dedoublet_adata.nf"
include { DEDOUBLET_ADATA as DEDOUBLET_RAW } from "../modules/dedoublet_adata.nf"
include { EXTRACT_EMBEDDING } from "../modules/extract_embedding.nf"


workflow DOUBLETS {
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

        DEDOUBLET_INTEGRATIONS(
            ch_integrations,
            SOLO.out
        )

        EXTRACT_EMBEDDING(DEDOUBLET_INTEGRATIONS.out)

        DEDOUBLET_RAW(
            ch_raw,
            SOLO.out
        )

    emit:
        solo = SOLO.out
        integrations = DEDOUBLET_INTEGRATIONS.out
        raw = DEDOUBLET_RAW.out
        obsm = EXTRACT_EMBEDDING.out
}