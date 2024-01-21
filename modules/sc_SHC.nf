process SCSHC_CLUSTERING {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.1"
    label "process_medium"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), val("${clustering_key}"), path("${meta.id}.${clustering_key}.clustering.h5ad")

    script:
    clustering_key = "scSHC"
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
    adata.obs["${clustering_key}"] = anndata2ri.rpy2py(result[0])
    adata.obs["${clustering_key}"] = adata.obs["${clustering_key}"].astype(int).astype(str).astype("category")

    adata.write_h5ad("${meta.id}.${clustering_key}.clustering.h5ad")
    """
}

process SCSHC_CLUSTERING_QC {
    tag "${meta.id}"

    container = "bigdatainbiomedicine/sc-rpy:1.1"
    publishDir "${params.outdir}/composition", mode: "${params.publish_mode}", pattern: "*.png"
    label "process_medium"

    input:
    tuple val(meta), val(clustering_key), path(adata)

    output:
    tuple val(meta), val(clustering_key), path("${meta.id}.${clustering_key}.clustering.qc.h5ad")

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
    clusters = anndata2ri.py2rpy(adata.obs['${clustering_key}'])

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

    adata.obs["qc"] = adata.obs["${clustering_key}"].map(value_mapping).fillna(0).astype(float)
    adata.write_h5ad("${meta.id}.${clustering_key}.clustering.qc.h5ad")

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

    DotExporter(root).to_picture("${meta.id}.${clustering_key}.qc.png")
    """
}