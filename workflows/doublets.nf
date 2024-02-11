include { PREPARE_SOLO } from "../modules/prepare_solo.nf"
include { SOLO } from "../modules/solo.nf"
include { DEDOUBLET_ADATA as DEDOUBLET_INTEGRATIONS } from "../modules/dedoublet_adata.nf"
include { DEDOUBLET_ADATA as DEDOUBLET_COUNTS } from "../modules/dedoublet_adata.nf"
include { EXTRACT_EMBEDDING } from "../modules/extract_embedding.nf"


workflow DOUBLETS {
    take:
        ch_adata
        ch_hvgs
        ch_scanvi_model
        ch_integrations
        ch_raw
    
    main:
        PREPARE_SOLO(
            ch_adata,
            ch_hvgs
        )

        ch_batches = PREPARE_SOLO.out.batches
            .splitText().flatten()
            .map{ batch -> batch.replace("\n", "") }

        SOLO(
            PREPARE_SOLO.out.adata.collect(),
            ch_scanvi_model.collect(),
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