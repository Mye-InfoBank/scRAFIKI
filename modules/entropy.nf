process ENTROPY {
  tag "${meta.id}"

  container "bigdatainbiomedicine/sc-rpy:1.0"
  label "process_medium"
  errorStrategy { task.attempt < 7 ? 'retry' : 'ignore'}

  input:
  tuple val(meta), path(adata)
  val(initial_smoothness)
  
  output:
  tuple val(meta), path("${meta.id}.entropy.pkl")

  when:
  task.ext.when == null || task.ext.when
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import anndata2ri
  import rpy2.robjects as ro
  import numpy as np
  rogue = ro.packages.importr('ROGUE')
  _ = ro.packages.importr('tibble')

  adata = ad.read_h5ad("$adata")
  sce = anndata2ri.py2rpy(adata)
  expression = ro.r("assay")(sce, "counts")

  label_col = "${meta.id}"
  sample_col = "batch"

  labels = anndata2ri.py2rpy(adata.obs[label_col].astype(str))
  samples = anndata2ri.py2rpy(adata.obs[sample_col])

  smoothness = ${initial_smoothness} * (1.2 ** (${task.attempt} - 1))

  result = rogue.rogue(expression, labels=labels, samples=samples, platform="UMI", span=smoothness)
  result = anndata2ri.rpy2py(result)
  entropy_dict = result.to_dict()

  df_entropy = adata.obs.apply(
    lambda row: entropy_dict.get(row[label_col], {}).get(row[sample_col], np.nan), axis=1
  )

  df_entropy.to_pickle("${meta.id}.entropy.pkl")
  """
}