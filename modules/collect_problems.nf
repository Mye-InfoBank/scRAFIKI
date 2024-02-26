process COLLECT_PROBLEMS {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  publishDir "${params.outdir}", mode: "${params.publish_mode}"

  label "process_single"

  input:
  path(problems)

  output:
  path("problems.txt")
  
  script:
  """
  #!/opt/conda/bin/python

  import os

  paths = "${problems.join(' ')}".split(' ')

  concatenated = ""
  suffix = ".problems.txt"
  for path in paths:
    name = os.path.basename(path)[0:-len(suffix)]
    concatenated += name + "\\n"
    with open(path, 'r') as f:
      lines = f.readlines()
    # Append, use tab indentation
    for line in lines:
      concatenated += "\\t" + line
    concatenated += "\\n"
  
  with open("problems.txt", 'w') as f:
    f.write(concatenated)
  """
}