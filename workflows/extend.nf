// Workflows
include { PREPROCESSING } from "../subworkflows/preprocessing.nf"


workflow EXTEND {
    if (!params.samplesheet) { 
        exit 1, 'Samplesheet not specified!'
    }

    ch_samplesheet = file(params.samplesheet)

    if (!params.base) {
        exit 1, 'Base path not specified!'
    }

    ch_base = Channel.value(file(params.base)).map{ base -> [[id: "base"], base]}

    if (!params.model) {
        exit 1, 'Model not specified!'
    }

    ch_model = Channel.value(file(params.model)).map{ model -> [[id: "model"], model]}

    PREPROCESSING(ch_samplesheet, ch_base)

    ch_adata = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obs = Channel.empty()
    
    emit:
        adata = ch_adata
        obsm = ch_obsm
        obs = ch_obs
}