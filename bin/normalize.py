#!/opt/conda/bin/python

#Python packages
import anndata as ad
from scipy.sparse import csr_matrix
import argparse

import anndata2ri
import rpy2.robjects as ro
seurat = ro.packages.importr('Seurat')

parser = argparse.ArgumentParser(description='Normalization of anndata object')
parser.add_argument('--input', type=str, help='input anndata object')
parser.add_argument('--output', type=str, help='output anndata object')
args = parser.parse_args()

adata = ad.read_h5ad(args.input)
sc_experiment = anndata2ri.py2rpy(adata)

# Convert to Seurat object using the as.Seurat function
# Need to extract the function beforehand, since it contains a dot
as_seurat = ro.r('as.Seurat')
seurat_object = as_seurat(sc_experiment, counts="ambient", data=ro.r('NULL'))

# transformed is a Seurat object with an assay called "SCT"
# https://github.com/satijalab/seurat/issues/4296
transformed = seurat.SCTransform(seurat_object, assay="originalexp", min_cells=0, method="glmGamPoi")

# Convert to singleCellExperiment object
as_sce = ro.r('as.SingleCellExperiment')
transformed_sce = as_sce(transformed)

# Convert to annData object
transformed_adata = anndata2ri.rpy2py(transformed_sce)

# Add the normalized matrix to the anndata object
adata.layers["normalized"] = csr_matrix(transformed_adata.X)
adata.X = adata.layers["normalized"]

# Write the anndata object to a file
adata.write_h5ad(args.output)