process BENCHMARK_INTEGRATIONS {
    tag "$meta2.id"

    label "process_high"

    container = "bigdatainbiomedicine/sc-scib:1.3.1"

    input: 
        tuple val(meta1), path(uncorrected)
        tuple val(meta2), path(integrated), val(integration_type)
        val(hvgs)

    output:
        path("${meta2.integration}.csv")
    
    when:
    task.ext.when == null || task.ext.when

    script:
        """
        scIB.py -u ${uncorrected} -i ${integrated} \
         -m ${meta2.integration} -o ${meta2.integration}.csv -b batch -l cell_type \
         --type ${integration_type} --hvgs ${hvgs}
        """
}

process MERGE_BENCHMARKS {
    label "process_single"

    publishDir "${params.outdir}", mode: "${params.publish_mode}"

    container = "bigdatainbiomedicine/sc-rpy:1.2"

    input:
        path(benchmarks)

    output:
        path("benchmarking.tsv")

    script:
        """
        #!/opt/conda/bin/python

        import pandas as pd

        benchmark_paths = "${benchmarks}".split(" ")
        benchmarks = [pd.read_csv(path, index_col=0) for path in benchmark_paths]

        merged = pd.concat(benchmarks, axis=1)

        merged.to_csv("benchmarking.tsv", sep="\\t")
        """
}