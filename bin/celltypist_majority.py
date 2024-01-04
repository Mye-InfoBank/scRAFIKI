#!/opt/conda/bin/python

import anndata as ad
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_clustering", type=str, help="input clustering h5ad file")
parser.add_argument("--input_celltypist", type=str, help="input celltypist h5ad file")
parser.add_argument("--cluster", type=str, help="cluster column name", default="leiden")
parser.add_argument("--min_prob", type=float, help="minimum probability", default=0)
parser.add_argument("--output", type=str, help="output h5ad file")

args = parser.parse_args()

adata = ad.read_h5ad(args.input_clustering)
adata_celltypist = ad.read_h5ad(args.input_celltypist)

predictions = adata_celltypist.obs["celltypist_prediction"]
clustering = adata.obs[args.cluster]

# The following section is based on https://github.com/Teichlab/celltypist/blob/777f8e92085492fcf61ff4a776c2a83600d7d1c0/celltypist/classifier.py#L444C10-L478C10

votes = pd.crosstab(predictions, clustering)
majority = votes.idxmax(axis=0).astype(str)

freqs = (votes / votes.sum(axis=0).values).max(axis=0)
majority[freqs < args.min_prob] = "Heterogeneous"
majority = majority[clustering].reset_index()
majority.index = predictions.index
majority.columns = ['clustering', 'majority_voting']
majority['majority_voting'] = majority['majority_voting'].astype('category')

# End of section

adata.obs["celltypist_majority"] = majority['majority_voting']

adata.write_h5ad(args.output)