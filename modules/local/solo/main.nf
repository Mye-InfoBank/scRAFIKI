#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process SOLO {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-python"

    label "process_medium"

    input:
        tuple val(meta), path(adata)
        tuple val(meta2), path(scvi_model)

    output:
        tuple val(meta), path("solo_${meta.id}.tsv")

    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    from threadpoolctl import threadpool_limits
    threadpool_limits(${task.cpus})

    import scvi

    def set_all_seeds(seed=0):
        import os
        import random
        import numpy as np
        import torch

        scvi.settings.seed = seed
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
        os.environ["PYTHONHASHSEED"] = str(seed)  # Python general
        np.random.seed(seed)  # Numpy random
        random.seed(seed)  # Python random

        torch.manual_seed(seed)
        torch.use_deterministic_algorithms(True)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)  # For multiGPU

    set_all_seeds()

    adata = sc.read_h5ad("${adata}")

    scvi.model.SCANVI.setup_anndata(adata, batch_key="batch", labels_key="cell_type", unlabeled_category="Unknown")
    scvi_model = scvi.model.SCANVI.load("${scvi_model}", adata=adata)
    solo = scvi.external.SOLO.from_scvi_model(scvi_model)
    solo.train()
    res = solo.predict()
    res["label"] = solo.predict(False)

    res.to_csv("solo_${meta.id}.tsv", sep="\\t")
    """
}

process FILTER_SOLO {
    container "bigdatainbiomedicine/sc-python"
    input:
        tuple val(batch), path(adata), path(info_df)

    output:
        tuple val(batch), path("*_dedup.h5ad")

    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    import pandas as pd
    from threadpoolctl import threadpool_limits
    threadpool_limits(${task.cpus})

    adata = sc.read_h5ad("${adata}")
    info_df = pd.read_csv("${info_df}", sep="\\t", index_col=0)

    adata_filtered = adata[info_df['label'] == "singlet"].copy()

    adata_filtered.write("solo_${batch}_dedup.h5ad")
    """
}