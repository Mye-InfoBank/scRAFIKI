process BENCHMARK_INTEGRATIONS {
    tag "$meta2.id"

    label "process_high"
    label "scale_resources"

    container = "bigdatainbiomedicine/sc-scib"
    cpus = 4
    memory = {50.GB * task.attempt}
    maxRetries = 4
    errorStrategy = 'retry'
    
    input: 
        tuple val(meta1), path(uncorrected)
        tuple val(meta2), path(integrated), val(integration_type)
        val(organism)
        val(hvgs)

    output:
        path("${meta2.integration}.csv")

    script:
        """
        scIB.py -u ${uncorrected} -i ${integrated} \
         -m ${meta2.integration} -o ${meta2.integration}.csv -b batch -l cell_type --organism ${organism} \
         --type ${integration_type} -f --hvgs ${hvgs}
        """
}