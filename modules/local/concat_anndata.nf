process CONCAT_ADATA {
  container = "bigdatainbiomedicine/sc-python"
  input:
    file anndatas
  
  output:
    file("concatenated.h5ad")
  
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