process BENCHMARK_INTEGRATIONS {
    tag "$meta.id"

    label "process_high"
    label "scale_resources"

    container = "bigdatainbiomedicine/sc-python"
    
    input: 
        path(uncorrected)
        tuple val(meta), path(integrated), val(integration_type), val(embed_key)
        val(organism)
        val(hvgs)

    output:
        path("${meta.integration}.csv")

    script:
        """
        scIB.py -u ${uncorrected} -i ${integrated} \
         -m ${meta.integration} -o ${meta.integration}.csv -b batch -l celltype --organism ${organism} \
         --type ${integration_type} --embed_key ${embed_key} -f --hvgs ${hvgs}
        """
}