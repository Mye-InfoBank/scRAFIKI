process DROPLET_CORRECT {
  tag "${meta.id}"
  container "bigdatainbiomedicine/sc-python"

  label "process_medium"

  input:
  tuple val(meta), path(adata)
  tuple val(meta2), path(solo)
  tuple val(meta3), path(decontX)
  
  output:
  tuple val(meta), path("droplet_corrected.h5ad")
  
  script:
  """
  #!/usr/bin/env python

  import anndata as ad
  import pandas as pd

  adata = ad.read_h5ad("$adata")
  solo_df = pd.read_csv("$solo", sep="\t", index_col=0)
  decontX_adata = ad.read_h5ad("$decontX")

  adata.layers["decontX"] = decontX_adata.X
  adata.obs["solo_doublet_score"] = solo_df["doublet"].values
  adata.obs["solo_singlet_score"] = solo_df["singlet"].values
  adata.obs["solo_label"] = solo_df["label"].values

  adata.write("droplet_corrected.h5ad")
  """
}