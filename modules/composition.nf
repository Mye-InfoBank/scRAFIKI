process COMPOSITION {
    container "bigdatainbiomedicine/sc-rpy:1.1"

    publishDir "${params.outdir}/composition", mode: "${params.publish_mode}"

    label "process_medium"

    input:
    tuple val(meta), path(adata)
    
    output:
    tuple val(meta), path("*.png")

    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import scanpy as sc
    import seaborn as sns
    import matplotlib.pyplot as plt
    import pandas as pd
    import upsetplot

    sns.set_theme(style="whitegrid")
    sns.set_context("paper")
    sns.set_palette("colorblind")
    sns.set(font_scale=1.5)

    adata = ad.read_h5ad("${adata}")

    def plot(data: pd.DataFrame, comparison: str):
        fractions = data.groupby(["dataset", comparison], observed=False).size().unstack()
        fractions = fractions.div(fractions.sum(axis=1), axis=0)

        fractions.plot.barh(stacked=True, title=f"Dataset composition ({comparison.capitalize()})")
        plt.legend(bbox_to_anchor=(1, 1.03), loc="upper left")
        plt.savefig(f"composition:{comparison}.png", bbox_inches="tight")
    
    for col in ["sex", "cell_type", "condition", "tissue"]:
        plot(adata.obs, col)

    # Split into multiple adatas, based on dataset
    datasets = adata.obs["dataset"].unique()
    dataset_genes = {}
    for dataset in datasets:
        adata_dataset = adata[adata.obs["dataset"] == dataset].copy()
        # Keep only genes with at least 1 count in at least 1 cell
        sc.pp.filter_genes(adata_dataset, min_cells=1)
        dataset_genes[dataset] = adata_dataset.var_names

    plot_data = upsetplot.from_contents(dataset_genes)

    upsetplot.plot(plot_data, sort_by="cardinality", show_counts=True, min_subset_size=10)
    plt.savefig("upset:genes.png")
    """
}