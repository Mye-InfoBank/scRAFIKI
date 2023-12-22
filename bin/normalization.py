#!/usr/bin/env python3

#Python packages
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse
import rpy2.robjects as ro
import anndata2ri
import sys

#R packages
ro.r('library(Seurat)')
ro.r('library(scater)')
anndata2ri.activate()

adata_path = sys.argv[1]
output_file = sys.argv[2]

adata = sc.read_h5ad(adata_path)

sc.pp.filter_genes(adata, min_cells=5)
    
if issparse(adata.X):
    if not adata.X.has_sorted_indices:
        adata.X.sort_indices()

for key in adata.layers:
    if issparse(adata.layers[key]):
        if not adata.layers[key].has_sorted_indices:
            adata.layers[key].sort_indices()

ro.globalenv['adata'] = adata

ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL)')

ro.r('res <- SCTransform(object=seurat_obj, assay = "originalexp", return.only.var.genes = FALSE, do.correct.umi = FALSE)')

norm_x = ro.r('res@assays$SCT@scale.data').T

adata_transformed = ad.AnnData(X=norm_x, obs=adata.obs)
# BASICALLY number of features differ: 6660 (original) vs 6645 (normalized)
# -> cannot concat the two anndata objects, but have 2 different ones now
# this can be due to filtering and used normalisation method (eg quality control...)

#adata_transformed.layers['normalized_data'] = adata_transformed.X

adata_transformed.write_h5ad(output_file)
