include { DECONTX } from "../modules/decontX.nf"
include { NORMALIZE } from "../modules/normalize.nf"

workflow COUNTS {
    take:
        ch_preprocessed
        normalization

    main:
        DECONTX(ch_preprocessed)

        NORMALIZE(DECONTX.out, normalization)

    emit:
        NORMALIZE.out
}