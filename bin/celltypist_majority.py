#!/opt/conda/bin/python

import anndata as ad
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_clustering", type=str, help="input clustering pickle file")
parser.add_argument("--input_celltypist", type=str, help="input celltypist pickle file")
parser.add_argument("--clustering_key", type=str, help="cluster column name", default="leiden")
parser.add_argument("--min_prob", type=float, help="minimum probability", default=0)
parser.add_argument("--output", type=str, help="output pickle file")

args = parser.parse_args()

clustering = pd.read_pickle(args.input_clustering)
df_celltypist = pd.read_pickle(args.input_celltypist)

predictions = df_celltypist["celltypist_prediction"]

majority_key = args.clustering_key + "_celltypist_majority"

# The following section is based on https://github.com/Teichlab/celltypist/blob/777f8e92085492fcf61ff4a776c2a83600d7d1c0/celltypist/classifier.py#L444C10-L478C10

votes = pd.crosstab(predictions, clustering)
majority = votes.idxmax(axis=0).astype(str)

freqs = (votes / votes.sum(axis=0).values).max(axis=0)
majority[freqs < args.min_prob] = "Heterogeneous"
majority = majority[clustering].reset_index()
majority.index = predictions.index
majority.columns = ['clustering', majority_key]
majority[majority_key] = majority[majority_key].astype('category')

# End of section

majority[majority_key].to_pickle(args.output)