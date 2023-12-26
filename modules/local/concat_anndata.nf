process CONCAT_ADATA {
  tag "$meta.id"
  container = "bigdatainbiomedicine/sc-python"
  label "process_medium"

  input:
    tuple val(meta), file(anndatas)
  
  output:
    tuple val(meta), file("concatenated.h5ad")
  
  script:
  """
    #!/usr/bin/env python
    import scanpy as sc
    import anndata as ad

    adata_list = [sc.read_h5ad(adata_path) for adata_path in ["${anndatas.join("\", \"")}"]]
    concatenated = ad.concat(adata_list)

    concatenated.write("concatenated.h5ad")
  """
}