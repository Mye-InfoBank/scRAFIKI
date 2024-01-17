process ENTROPY {
  tag "${meta.id}"

  container "bigdatainbiomedicine/sc-rpy"
  label "process_medium"

  input:
  tuple val(meta), val(clustering_key), path(adata)
  
  output:
  tuple val(meta), val(clustering_key), path("${meta.id}.${clustering_key}.entropy.h5ad")
  
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

  label_col = "${clustering_key}"
  sample_col = "batch"

  labels = anndata2ri.py2rpy(adata.obs[label_col].astype(str))
  samples = anndata2ri.py2rpy(adata.obs[sample_col])

  smoothness = 0.5 * (1.2 ** (${task.attempt} - 1))

  result = rogue.rogue(expression, labels=labels, samples=samples, platform="UMI", span=smoothness)
  result = anndata2ri.rpy2py(result)
  entropy_dict = result.to_dict()

  adata.obs["entropy"] = adata.obs.apply(
    lambda row: entropy_dict.get(row[label_col], {}).get(row[sample_col], np.nan), axis=1
  )

  adata.write_h5ad("${meta.id}.${clustering_key}.entropy.h5ad")
  """
}