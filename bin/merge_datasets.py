#!/opt/conda/bin/python

import argparse
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
import numpy as np

columns_required = {
    "sex": False,
    "batch": True,
    "cell_type": False,
    "condition": False,
    "patient": True,
    "tissue": True,
    "dataset": True,
    "transfer": True
}

parser = argparse.ArgumentParser(description="Merge datasets")
parser.add_argument("--input", help="Input file", type=str, nargs="+")
parser.add_argument("--output_integration", help="Output file containing only cells which do not require transfer learning", type=str)
parser.add_argument("--output_intersection", help="Output file containing all cells but gene intersection", type=str)
parser.add_argument("--output_transfer", help="Output file containing all cells which require transfer learning", type=str)
parser.add_argument("--output_counts", help="Output file, outer join of cells and genes", type=str)
parser.add_argument("--min_cells", help='Minimum number of cells to keep a gene', type=int, required=False, default=50)
parser.add_argument("--custom_metadata", help="Additional metadata columns to include", type=str, nargs="*")
parser.add_argument("--custom_genes", help="Additional genes to include", type=str, nargs="*")

args = parser.parse_args()

columns_additional = {column: False for column in args.custom_metadata}
columns_considered = {**columns_required, **columns_additional}

datasets = [ad.read_h5ad(f) for f in args.input]

for file_name, dataset in zip(args.input, datasets):
    # Make sure dataset is not empty
    if dataset.shape[0] == 0 or dataset.shape[1] == 0:
        raise ValueError(
            f'Dataset is empty: {file_name}'
        )
    # Make sure all required columns are present
    for column, required in columns_considered.items():
        if column not in dataset.obs.columns:
            if required:
                raise ValueError(
                    f"Column {column} is required but not found in {file_name}"
                )
            else:
                dataset.obs[column] = "Unknown"

    # Subset columns
    dataset.obs = dataset.obs[columns_considered.keys()]

adata = ad.concat(datasets)
adata_outer = ad.concat(datasets, join='outer')

additional_genes = [gene for gene in args.custom_genes if gene not in adata.var_names]

# Add custom genes from outer join to the intersection
if additional_genes:
    adata_additional = adata_outer[adata.obs_names, additional_genes]
    adata_concatenated = ad.concat([adata, adata_additional], join="outer", axis=1)
    adata_concatenated.obs, adata_concatenated.obsm = adata.obs, adata.obsm
    adata = adata_concatenated

# Convert to CSR matrix
adata.X = csr_matrix(adata.X)
adata_outer.X = csr_matrix(adata_outer.X)

# Make sure that there are no underscores in the cell names
adata.obs_names = adata.obs_names.str.replace("_", "-")
adata.obs_names_make_unique()
adata_outer.obs_names = adata.obs_names

# Filter genes with no counts in core atlas
gene_mask, _ = sc.pp.filter_genes(adata[~adata.obs["transfer"]], min_cells=1, inplace=False)
adata = adata[:, gene_mask]

# Filter cells with no counts
cell_mask, _ = sc.pp.filter_cells(adata, min_genes=1, inplace=False)
adata = adata[cell_mask, :]
adata_outer = adata_outer[cell_mask, :]

# Filter genes with too few occurrences in outer join
sc.pp.filter_genes(adata_outer, min_cells=args.min_cells)

adata.obs["batch"] = adata.obs["dataset"].astype(str) + "_" + adata.obs["batch"].astype(str)
adata.obs["patient"] = adata.obs["dataset"].astype(str) + "_" + adata.obs["patient"].astype(str)

def to_Florent_case(s: str):
    corrected = s.lower().strip()

    if corrected in ["na", "nan", "null", "unknown"]:
        return "Unknown"

    corrected = s \
        .replace(" ", "_") \
        .replace("-", "_")

    corrected = "".join([c if c.isalnum() or c == "_" else "" for c in corrected])

    # Make sure there is never more than one underscore
    corrected = corrected.replace("__", "_")

    if corrected.endswith("s"):
        corrected = corrected[:-1]

    corrected = corrected.strip(" _")

    if not corrected:
        return "Unknown"

    return corrected[0].upper() + corrected[1:]

for column in columns_considered.keys():
    if column == "transfer":
        continue
    # Convert first to string and then to category
    adata.obs[column] = adata.obs[column].astype(str).fillna("Unknown").apply(to_Florent_case).astype("category")

adata_outer.obs = adata.obs

adata.layers["counts"] = adata.X
adata_outer.layers["counts"] = adata_outer.X

if any(adata.obs["transfer"]):
    adata_transfer = adata[adata.obs["transfer"]]
    adata_transfer.write_h5ad(args.output_transfer)

adata_notransfer = adata[~adata.obs["transfer"]]
adata_notransfer.write_h5ad(args.output_integration)

adata.write_h5ad(args.output_intersection)
adata_outer.write_h5ad(args.output_counts)