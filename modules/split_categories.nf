process SPLIT_CATEGORIES {
    tag "$meta.id"

    container "bigdatainbiomedicine/sc-rpy:1.0"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    path(df_categories)
    val(category)
    
    output:
    tuple val(meta), path("*.category.h5ad"), emit: adata
    
    script:
    if (df_categories)
        """
        #!/opt/conda/bin/python

        import anndata as ad
        import pandas as pd

        adata = ad.read_h5ad("${adata}")
        df_categories = pd.read_csv("${df_categories}", index_col=0, comment="#")
        adata.obs["category"] = df_categories.loc[adata.obs.index, "${category}"]

        for category in adata.obs["category"].unique():
            adata_category = adata[adata.obs["category"] == category]
            adata_category.write_h5ad(f"{category.replace(' ', '_')}.category.h5ad")
        """
    else
        """
        #!/opt/conda/bin/python
        
        import anndata as ad

        adata = ad.read_h5ad("${adata}")

        for category in adata.obs["${category}"].unique():
            adata_category = adata[adata.obs["${category}"] == category]
            adata_category.write_h5ad(f"{category.replace(' ', '_')}.category.h5ad")
        """
}