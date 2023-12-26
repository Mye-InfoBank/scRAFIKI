process ADATA_METRICS {
    tag "${meta.id}"
    container "bigdatainbiomedicine/sc-python"

    label "process_medium"

    input:
        tuple val(meta), file(adata)
    
    output:
        file("${meta.id}_metrics.csv")
    
    script:
    """
        #!/usr/bin/env python3

        import scanpy as sc
        import pandas as pd

        adata = sc.read_h5ad("${adata}")

        patient_count = len(adata.obs['patient'].unique())
        cell_count = len(adata.obs_names)
        celltype_counts = adata.obs['celltype'].value_counts().to_dict() if 'celltype' in adata.obs.columns else {}

        celltypes = list(celltype_counts.keys())
        celltypes.sort()

        celltype_count_values = [celltype_counts[celltype] for celltype in celltypes]

        # Create a dataframe with columns: patient_count, cell_count, celltype_1, celltype_2, ...
        df = pd.DataFrame([[patient_count, cell_count] + celltype_count_values], columns=['patient_count', 'cell_count'] + celltypes)
        df.index = ["${meta.id}"]

        df.to_csv("${meta.id}_metrics.csv", index=True)
    """
}

process COMBINE_ADATA_METRICS {
    container = "bigdatainbiomedicine/sc-python"

    input:
        file(metrics)
    
    output:
        file("metrics.csv")
    
    script:
    """
        #!/usr/bin/env python

        import pandas as pd

        paths = "${metrics}".split()
        dataframes = [pd.read_csv(path, index_col=0) for path in paths]
        df = pd.concat(dataframes)

        df.to_csv("metrics.csv", index=True)
    """
}