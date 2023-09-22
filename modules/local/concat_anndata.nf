process CONCAT_ADATA {
  input:
    val anndatas
  
  output:
    file("concatenated.h5ad")
  
  script:
  """
    import scanpy as sc
    import anndata as ad

    anndata_list = [sc.read_h5ad(adata) for adata in ["${anndatas.join("\", \"")}}"]]
    concatenated = ad.concatenate(anndata_list)

    concatenated.write("concatenated.h5ad")
  """
}