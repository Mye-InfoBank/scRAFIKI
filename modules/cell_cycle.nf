process CELL_CYCLE {
    tag "$meta.id"

    label "process_medium"

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

        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        sc.pp.scale(adata)

        scib.pp.score_cell_cycle(adata, "${organism}")

        adata.write_h5ad("${meta.id}.cell_cycle.h5ad")
        """
}