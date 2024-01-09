include { DECONTX } from "../modules/decontX.nf"
include { NORMALIZE } from "../modules/normalize.nf"

workflow COUNTS {
    take:
        ch_preprocessed
        normalization

    main:
        DECONTX(ch_preprocessed)

        ch_counts = DECONTX.out

        if (params.normalize) {
            NORMALIZE(ch_counts, normalization)

            ch_counts = NORMALIZE.out
        }

    emit:
        ch_counts
}