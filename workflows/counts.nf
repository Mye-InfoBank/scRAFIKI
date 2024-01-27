include { DECONTX } from "../modules/decontX.nf"
include { NORMALIZE } from "../modules/normalize.nf"

workflow COUNTS {
    take:
        ch_preprocessed
        normalization

    main:
        if (params.decontX) {
            DECONTX(ch_preprocessed)
            ch_counts = DECONTX.out
        } else {
            ch_counts = ch_preprocessed
        }

        NORMALIZE(ch_counts, normalization)

    emit:
        NORMALIZE.out
}