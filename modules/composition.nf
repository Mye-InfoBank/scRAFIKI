process COMPOSITION {
    container "bigdatainbiomedicine/sc-rpy:1.0"

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
    """
}