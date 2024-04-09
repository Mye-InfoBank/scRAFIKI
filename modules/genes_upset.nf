process GENES_UPSET {
    container "bigdatainbiomedicine/sc-rpy:1.2"

    publishDir "${params.outdir}", mode: "${params.publish_mode}"

    label "process_medium"

    input:
    path(adatas)
    
    output:
    path("*.png"), optional: true

    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import scanpy as sc
    from collections import defaultdict
    import matplotlib.pyplot as plt
    import upsetplot

    adata_paths = "${adatas}".split(" ")
    dataset_genes = defaultdict(set)

    for adata_path in adata_paths:
        adata = ad.read_h5ad(adata_path)
        # Split into multiple adatas, based on dataset
        datasets = adata.obs["dataset"].unique()
        for dataset in datasets:
            adata_dataset = adata[adata.obs["dataset"] == dataset].copy()
            # Keep only genes with at least 1 count in at least 1 cell
            sc.pp.filter_genes(adata_dataset, min_cells=1)
            dataset_genes[dataset].update(adata_dataset.var_names)

    if not len(dataset_genes) > 1:
        exit()

    plot_data = upsetplot.from_contents(dataset_genes)

    upsetplot.plot(plot_data, sort_by="cardinality", show_counts=True, min_subset_size=10)
    plt.savefig("upset:genes.png")
    """
}