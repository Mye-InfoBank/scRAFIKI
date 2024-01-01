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
        cell_type_counts = adata.obs['cell_type'].value_counts().to_dict() if 'cell_type' in adata.obs.columns else {}

        cell_types = list(cell_type_counts.keys())
        cell_types.sort()

        cell_type_count_values = [cell_type_counts[cell_type] for cell_type in cell_types]

        # Create a dataframe with columns: patient_count, cell_count, cell_type_1, cell_type_2, ...
        df = pd.DataFrame([[patient_count, cell_count] + cell_type_count_values], columns=['patient_count', 'cell_count'] + cell_types)
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