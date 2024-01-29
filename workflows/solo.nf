include { SOLO } from "../modules/solo.nf"

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

    emit:
        SOLO.out
}