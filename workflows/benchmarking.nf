include { BENCHMARK_INTEGRATIONS } from "../modules/scIB.nf"
include { MERGE_BENCHMARKS } from "../modules/scIB.nf"

workflow BENCHMARKING {
    take:
        ch_uncorrected
        ch_integration
        hvgs
    
    main:
        BENCHMARK_INTEGRATIONS(ch_uncorrected, ch_integration, hvgs)
        MERGE_BENCHMARKS(BENCHMARK_INTEGRATIONS.out.collect())
}