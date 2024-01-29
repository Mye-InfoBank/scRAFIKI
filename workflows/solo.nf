include { SOLO } from "../modules/solo.nf"
include { MERGE_SOLO } from "../modules/solo.nf"

workflow W_SOLO {
    take:
        ch_hvgs
        ch_scanvi_model
        ch_batches
    
    main:
        SOLO(
            ch_hvgs,
            ch_scanvi_model,
            ch_batches
        )

        MERGE_SOLO(
            ch_hvgs,
            SOLO.out.map{ meta, adata -> adata}.collect()
        )

    emit:
        MERGE_SOLO.out
}