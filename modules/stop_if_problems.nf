process STOP_IF_PROBLEMS {
  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_single"

  errorStrategy 'terminate'

  input:
  path(problems)
  
  script:
  """
  echo "Problems found in input datasets. View problems.txt in the output directory for more details."
  exit 1
  """
}