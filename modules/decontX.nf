process DECONTX {
  tag "$meta.id"

  container "bigdatainbiomedicine/sc-rpy:1.0"

  label "process_medium"
  label "process_high_memory"

  input:
    tuple val(meta), file(adata)
  
  output:
    tuple val(meta), file("${meta.id}.ambient.h5ad")
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import anndata2ri
  import rpy2.robjects as ro
  celda = ro.packages.importr('celda')

  adata = ad.read_h5ad("${adata}")
  sc_experiment = anndata2ri.py2rpy(adata)

  corrected = celda.decontX(sc_experiment, batch=adata.obs['batch'].tolist())
  counts = celda.decontXcounts(corrected)

  adata.layers['ambient'] = anndata2ri.rpy2py(counts).T
  adata.X = adata.layers['ambient']
  adata.write_h5ad("${meta.id}.ambient.h5ad")
  """
}