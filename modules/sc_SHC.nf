process SCSHC_CLUSTERING {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.2"
    label "process_high"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(new_meta), path("${new_meta.id}.clustering.h5ad"), emit: adata
    tuple val(new_meta), path("${new_meta.id}.clustering.pkl"), emit: table

    when:
    task.ext.when == null || task.ext.when

    script:
    new_meta = meta + [id: "${meta.integration}_scSHC"]
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import anndata2ri
    import rpy2
    import rpy2.robjects as ro
    scSHC = ro.packages.importr('scSHC')

    adata = ad.read_h5ad("${adata}")
    sce = anndata2ri.py2rpy(adata)
    expression = ro.r("counts")(sce)
    batches = anndata2ri.py2rpy(adata.obs['batch'])

    n_genes = min(adata.shape[1], 2500)

    result = scSHC.scSHC(expression, batch=batches, cores=${task.cpus}, num_features=n_genes)
    adata.obs["${new_meta.id}"] = anndata2ri.rpy2py(result[0])
    adata.obs["${new_meta.id}"] = adata.obs["${new_meta.id}"].astype(int).astype(str).astype("category")

    adata.obs["${new_meta.id}"].to_pickle("${new_meta.id}.clustering.pkl")
    adata.write_h5ad("${new_meta.id}.clustering.h5ad")
    """
}

process SCSHC_CLUSTERING_QC {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.2"
    publishDir "${params.outdir}/clustering_qc", mode: "${params.publish_mode}", pattern: "*.png"
    label "process_high"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("${meta.id}.qc.pkl"), emit: table
    tuple val(meta), path("${meta.id}.qc.png"), emit: png, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/opt/conda/bin/python

    import anndata as ad
    import anndata2ri
    from anytree import AnyNode
    from anytree.exporter import DotExporter
    import pandas as pd
    import rpy2
    import rpy2.robjects as ro
    scSHC = ro.packages.importr('scSHC')
    datatree = ro.packages.importr('data.tree')

    adata = ad.read_h5ad("${adata}")
    sce = anndata2ri.py2rpy(adata)
    expression = ro.r("counts")(sce)
    batches = anndata2ri.py2rpy(adata.obs['batch'])
    clusters = anndata2ri.py2rpy(adata.obs['${meta.id}'])

    n_genes = min(adata.shape[1], 2500)

    result = scSHC.testClusters(expression, clusters, batch=batches, cores=${task.cpus}, num_features=n_genes)

    network_df = anndata2ri.rpy2py(datatree.ToDataFrameNetwork(result[1]))
    all_nodes = set(network_df["from"].to_list() + network_df["to"].to_list())

    name_mapping = {
        node: node.split(":")[0] for node in all_nodes
    }

    splitted_nodes = [node.split(":") for node in all_nodes]
    value_mapping = {
        splitted_node[0].split(" ")[1] : splitted_node[1]
        for splitted_node in splitted_nodes if len(splitted_node) == 2 and splitted_node[0].startswith("Cluster")
    }

    df_qc = adata.obs["${meta.id}"].map(value_mapping).fillna(0).astype(float)
    # Write to tsv
    df_qc.to_pickle("${meta.id}.qc.pkl")

    if len(network_df) == 0:
        print("No clusters found")
        exit(0)

    root = AnyNode(id=network_df.iloc[0]["from"], name=network_df.iloc[0]["from"])
    nodes = {
        network_df.iloc[0]["from"]: root
    }

    for index, row in network_df.iterrows():
        if row["from"] not in nodes:
            raise Exception("Node not found")
        if row["to"] not in nodes:
            nodes[row["to"]] = AnyNode(id=row["to"], name=row["to"], parent=nodes[row["from"]])

    DotExporter(root).to_picture("${meta.id}.qc.png")
    """
}