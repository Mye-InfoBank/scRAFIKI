include { SPLIT_BATCHES } from "../modules/split_batches.nf"
include { DECONTX } from "../modules/decontX.nf"
include { CONCAT_ADATA as CONCAT_DECONTX } from "../modules/concat_anndata.nf"
include { NORMALIZE } from "../modules/normalize.nf"

workflow COUNTS {
    take:
        ch_preprocessed

    main:
        SPLIT_BATCHES(ch_preprocessed)

        ch_batches = SPLIT_BATCHES.out
            .transpose()
            .map{ meta, adata -> 
                batch = adata.simpleName;
                return [[id: batch], adata]
            }
        
        DECONTX(ch_batches)

        CONCAT_DECONTX(
            DECONTX.out.map{ it[1] }.collect().map{ [[id: "counts"], it] }
        )

        NORMALIZE(CONCAT_DECONTX.out)

    emit:
        NORMALIZE.out
}