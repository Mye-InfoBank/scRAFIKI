process CELL_CYCLE {
    tag "$meta.id"

    label "process_high"
    label "process_high_memory"
    errorStrategy 'retry'

    container = "bigdatainbiomedicine/sc-scib:1.3"

    input:
        tuple val(meta), path(adata)
        val(organism)

    output:
        tuple val(meta), path("${meta.id}.cell_cycle.pkl")

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        #!/usr/bin/env python

        import scanpy as sc
        import scib

        print("Reading data")
        adata = sc.read_h5ad("${adata}")

        print("Preprocessing")
        print("\\tNormalizing")
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        print("\\tLogarithmizing")
        sc.pp.log1p(adata)
        print("\\tScaling")
        sc.pp.scale(adata)

        print("Cell cycle scoring")
        scib.pp.score_cell_cycle(adata, "${organism}")

        df_cell_cycle = adata.obs[["G2M_score", "S_score", "phase"]]
        df_cell_cycle.columns = ["G2M_score", "S_score", "cycle_phase"]

        print("Saving results")
        df_cell_cycle.to_pickle("${meta.id}.cell_cycle.pkl")
        """
}