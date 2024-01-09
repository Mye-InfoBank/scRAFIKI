#!/opt/conda/bin/python

import argparse
import anndata as ad
import scanpy as sc
from scipy.sparse import csc_matrix

columns_required = {
    "sex": False,
    "batch": True,
    "cell_type": False,
    "condition": False,
    "patient": True,
    "tissue": True,
    "dataset": True,
}

parser = argparse.ArgumentParser(description="Merge datasets")
parser.add_argument("--input", help="Input file", type=str, nargs="+")
parser.add_argument("--output", help="Output file", type=str)

args = parser.parse_args()

datasets = [ad.read_h5ad(f) for f in args.input]

for dataset in datasets:
    # Make sure all required columns are present
    for column, required in columns_required.items():
        if column not in dataset.obs.columns:
            if required:
                raise ValueError(
                    f"Column {column} is required but not found in {dataset}"
                )
            else:
                dataset.obs[column] = "unknown"

    # Subset columns
    dataset.obs = dataset.obs[columns_required.keys()]

merged = ad.concat(datasets, join="outer")

# Make sure that there are no underscores in the cell names
merged.obs_names = merged.obs_names.str.replace("_", "-")
merged.obs_names_make_unique()

# Make sure that there are no underscores in the gene names
merged.var_names = merged.var_names.str.replace("_", "-")
merged.var_names_make_unique()

# Perform minimal filtering to prevent NaNs
sc.pp.filter_cells(merged, min_genes=1)
sc.pp.filter_genes(merged, min_cells=1)

# Convert to CSC matrix
merged.X = csc_matrix(merged.X)

merged.obs["batch"] = merged.obs["dataset"].astype(str) + "_" + merged.obs["batch"].astype(str)
merged.obs["patient"] = merged.obs["dataset"].astype(str) + "_" + merged.obs["patient"].astype(str)

for column in columns_required.keys():
    # Convert first to string and then to category
    merged.obs[column] = merged.obs[column].astype(str).astype("category")

merged.layers["counts"] = merged.X.copy()

merged.write_h5ad(args.output)