process INTEGRATE {
  tag "${method}"
  container "bigdatainbiomedicine/sc-scib:1.2"

  label "process_high"

  input:
  tuple val(meta), path(input)
  tuple val(meta2), path(hvgs)
  val(method)
  
  output:
  tuple val(meta_out), path("${method}.h5ad"), emit: integrated
  tuple val(meta_out), path("model"), emit: model, optional: true

  script:
  meta_out = ["id": "${method}", "integration": "${method}"]
  """
  integrate.py --input ${input} --hvgs ${hvgs} --method ${method} --output ${method}.h5ad --cpus ${task.cpus}
  """
}