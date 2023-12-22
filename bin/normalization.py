#!/usr/bin/env python3

#Python packages
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
import rpy2.robjects as ro
import anndata2ri
import sys
import os
os.environ['NUMBA_CACHE_DIR'] = '/tmp/numba_cache'

#R packages
ro.r('library(Seurat)')
ro.r('library(scater)')
anndata2ri.activate()

adata_path = sys.argv[1]
output_file = sys.argv[2]

adata = sc.read_h5ad(adata_path)

gene_counts = adata.var_names.value_counts()
duplicates = gene_counts[gene_counts > 1]
if not duplicates.empty:
    adata.var_names_make_unique()

sc.pp.filter_genes(adata, min_cells=5)
    
if issparse(adata.X):
    if not adata.X.has_sorted_indices:
        adata.X.sort_indices()

for key in adata.layers:
    if issparse(adata.layers[key]):
        if not adata.layers[key].has_sorted_indices:
            adata.layers[key].sort_indices()

adata.obs_names_make_unique()

ro.globalenv['adata'] = adata

ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL)')

ro.r('res <- SCTransform(object=seurat_obj, assay = "originalexp", return.only.var.genes = FALSE, do.correct.umi = FALSE)')

norm_x = ro.r('res@assays$SCT@scale.data').T
norm_x = csr_matrix(norm_x)

adata_transformed = ad.AnnData(X=norm_x, obs=adata.obs)
# BASICALLY number of features differ: 6660 (original) vs 6645 (normalized)
# -> cannot concat the two anndata objects, but have 2 different ones now
# this can be due to filtering and used normalisation method (eg quality control...)

#adata_transformed.layers['normalized_data'] = adata_transformed.X

adata_transformed.write_h5ad(output_file)
