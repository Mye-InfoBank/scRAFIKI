process BENCHMARK_INTEGRATIONS {
    container = "bigdatainbiomedicine/sc-python"
    cpus = 4
    memory = {50.GB * task.attempt}
    maxRetries = 4
    errorStrategy = 'retry'
    
    input: 
        tuple path(uncorrected), val(name), val(embed_acc), path(integrated), val(integration_type)
        val(organism)

    output:
        path("${name}.csv")

    script:
        """
        scIB.py -u ${uncorrected} -i ${integrated} \
         -m ${name} -o ${name}.csv -b batch -l cell_type --organism ${organism} \
         --type ${integration_type} --embed_key ${embed_acc}
        """
}