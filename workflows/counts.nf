include { DECONTX } from "../modules/decontX.nf"
include { CELLBENDER } from "../modules/cellbender.nf"
include { H5_TO_H5AD } from "../modules/h5_to_h5ad.nf"
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
            CELLBENDER(ch_preprocessed)
            H5_TO_H5AD(CELLBENDER.out)
            ch_counts = H5_TO_H5AD.out
        }

        NORMALIZE(ch_counts, normalization)

    emit:
        NORMALIZE.out
}