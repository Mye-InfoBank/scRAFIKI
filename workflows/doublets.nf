include { SOLO } from "../modules/solo.nf"
include { DEDOUBLET_ADATA as DEDOUBLET_INTEGRATIONS } from "../modules/dedoublet_adata.nf"
include { DEDOUBLET_ADATA as DEDOUBLET_COUNTS } from "../modules/dedoublet_adata.nf"
include { EXTRACT_EMBEDDING } from "../modules/extract_embedding.nf"


workflow DOUBLETS {
    take:
        ch_hvgs
        ch_scanvi_model
        ch_integrations
        ch_raw
    
    main:
        SOLO(
            ch_hvgs,
            ch_scanvi_model
        )

        DEDOUBLET_INTEGRATIONS(
            ch_integrations,
            SOLO.out
        )

        EXTRACT_EMBEDDING(DEDOUBLET_INTEGRATIONS.out)

        DEDOUBLET_COUNTS(
            ch_raw,
            SOLO.out
        )

    emit:
        solo = SOLO.out
        integrations = DEDOUBLET_INTEGRATIONS.out
        counts = DEDOUBLET_COUNTS.out
        obsm = EXTRACT_EMBEDDING.out
}