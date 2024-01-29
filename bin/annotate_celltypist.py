#!/opt/conda/bin/python

import scanpy as sc
import warnings
import celltypist
from celltypist import models
import argparse

warnings.filterwarnings("ignore", category=FutureWarning)

parser = argparse.ArgumentParser()

parser.add_argument("--input", type=argparse.FileType("r"), help="input h5ad file")
parser.add_argument("--output", type=argparse.FileType("w"), help="output pickle file")
parser.add_argument("--model", type=str, help="The celltypist model file to use")

args = parser.parse_args()


adata = sc.read_h5ad(args.input.name)

adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(
    adata_celltypist, counts_per_cell_after=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()

models.download_models(model=args.model)
model = models.Model.load(args.model)
 
predictions = celltypist.annotate(
    adata_celltypist, model=model
)
predictions_adata = predictions.to_adata()

df_celltypist = predictions_adata.obs.loc[
    adata.obs.index, ["predicted_labels", "conf_score"]
]

df_celltypist.columns = ["cell_type:celltypist", "celltypist_confidence"]

df_celltypist.to_pickle(args.output.name)
