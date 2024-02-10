process INTEGRATE_SCANVI {
  tag "${method}"
  container "bigdatainbiomedicine/sc-scib:1.2"

  label "process_high"

  publishDir "${params.outdir}", mode: "${params.publish_mode}", pattern: "model"

  input:
  tuple val(meta), path(input)
  tuple val(meta2), path(scvi_model, stageAs: "scvi_model")
  
  output:
  tuple val(meta_out), path("${method}.h5ad"), emit: integrated
  tuple val(meta_out), path("model"), emit: model
  tuple val(meta_out), path("scanvi_labels.pkl"), emit: labels
  
  script:
  method = "scanvi"
  meta_out = ["id": "${method}", "integration": "${method}"]
  """
  integrate.py --input ${input} --method ${method} --scvi_model ${scvi_model} --output ${method}.h5ad
  """
}