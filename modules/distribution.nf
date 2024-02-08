process DISTRIBUTION {
    container "bigdatainbiomedicine/sc-rpy:1.0"

    publishDir "${params.outdir}/distribution", mode: "${params.publish_mode}"

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
        fractions = data.groupby(["cell_type", comparison], observed=False).size().unstack()
        fractions = fractions.div(fractions.sum(axis=1), axis=0)

        fractions.plot.barh(stacked=True, title=f"Cell type distribution ({comparison.capitalize()})")
        plt.legend(bbox_to_anchor=(1, 1.03), loc="upper left")
        plt.savefig(f"distribution:{comparison}.png", bbox_inches="tight")
    
    for col in ["sex", "batch", "condition", "tissue", "dataset"]:
        plot(adata.obs, col)
    """
}