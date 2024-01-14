process CELL_CYCLE {
    tag "$meta.id"

    label "process_high"
    label "process_high_memory"
    errorStrategy 'retry'

    container = "bigdatainbiomedicine/sc-scib"

    input:
        tuple val(meta), path(adata)
        val(organism)

    output:
        tuple val(meta), path("${meta.id}.cell_cycle.h5ad")

    script:
        """
        #!/usr/bin/env python

        import scanpy as sc
        import scib

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

        adata.write_h5ad("${meta.id}.cell_cycle.h5ad")
        """
}