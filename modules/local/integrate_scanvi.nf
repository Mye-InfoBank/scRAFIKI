process INTEGRATE_SCANVI {
  tag "${method}"
  container "bigdatainbiomedicine/sc-python"

  label "process_high"
  label "scale_resources"

  input:
  tuple val(meta), path(input)
  tuple val(meta2), path(scvi_model, stageAs: "scvi_model")
  
  output:
  tuple val(meta_out), path("${method}.h5ad"), emit: integrated
  tuple val(meta_out), path("model"), emit: model
  
  script:
  method = "scanvi"
  meta_out = ["id": "${method}", "integration": "${method}"]
  """
  integrate.py --input ${input} --method ${method} --scvi_model ${scvi_model} --output ${method}.h5ad
  """
}