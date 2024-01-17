process MERGE {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-rpy:1.0"

  publishDir "${params.outdir}", mode: "${params.publish_mode}"

  label "process_high_memory"

  input:
    tuple val(meta), path(original_adata, stageAs: 'input.h5ad')
    val (integration_names)
    path(integration_adatas)
    tuple val(meta2), path(solo)
    tuple val(meta3), path(counts)
    tuple val(meta4), path(celltypist)
    tuple val(meta5), val(cell_cycle)
    val(resolutions)
  
  output:
    file "merged.h5ad"
  
  script:
  """
  #!/opt/conda/bin/python

  import anndata as ad
  import pandas as pd

  adata = ad.read_h5ad("${original_adata}")

  for integration_name, integration_adata_path in zip(
        ["${integration_names.join("\",\"")}"],
        ["${integration_adatas.join("\",\"")}"]):
      integration_adata = ad.read_h5ad(integration_adata_path)
      adata.obsm['emb_' + integration_name] = integration_adata.obsm["X_emb"].copy()
      adata.obsm['X_' + integration_name] = integration_adata.obsm['X_umap'].copy()

      for resolution in ${resolutions}:
        res = float(resolution)
        leiden_key = f"leiden_{res:.2f}"
        majority_key = f"{leiden_key}_celltypist_majority"
        entropy_key = f"{leiden_key}_entropy"
        adata.obs[integration_name + '_' + leiden_key] = integration_adata.obs[leiden_key].copy()
        adata.obs[integration_name + '_' + majority_key] = integration_adata.obs[majority_key].copy()

        if entropy_key in integration_adata.obs.keys():
          adata.obs[integration_name + '_' + entropy_key] = integration_adata.obs[entropy_key].copy()

      del integration_adata

  solo_df = pd.read_csv("$solo", sep="\\t", index_col=0)
  counts_adata = ad.read_h5ad("$counts")

  for layer in counts_adata.layers.keys():
    adata.layers[layer] = counts_adata.layers[layer]
  adata.X = counts_adata.X
  del counts_adata

  celltypist_adata = ad.read_h5ad("$celltypist")
  adata.obs["celltypist_prediction"] = celltypist_adata.obs["celltypist_prediction"].values
  adata.obs["celltypist_conf_score"] = celltypist_adata.obs["celltypist_conf_score"].values
  del celltypist_adata

  cell_cycle_adata = ad.read_h5ad("$cell_cycle")
  adata.obs["S_score"] = cell_cycle_adata.obs["S_score"].values
  adata.obs["G2M_score"] = cell_cycle_adata.obs["G2M_score"].values
  adata.obs["cycle_phase"] = cell_cycle_adata.obs["phase"].values
  del cell_cycle_adata

  adata.obs["solo_doublet_score"] = solo_df["doublet"].values
  adata.obs["solo_singlet_score"] = solo_df["singlet"].values
  adata.obs["solo_label"] = solo_df["label"].values

  adata.write('merged.h5ad')
  """
}