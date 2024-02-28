#!/opt/conda/bin/python

import argparse
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
import numpy as np



parser = argparse.ArgumentParser(description="Merge datasets")
parser.add_argument("--input", help="Input file", type=str, nargs="+")
parser.add_argument("--output_integration", help="Output file containing only cells which do not require transfer learning", type=str)
parser.add_argument("--output_intersection", help="Output file containing all cells but gene intersection", type=str)
parser.add_argument("--output_transfer", help="Output file containing all cells which require transfer learning", type=str)
parser.add_argument("--output_counts", help="Output file, outer join of cells and genes", type=str)
parser.add_argument("--min_cells", help='Minimum number of cells to keep a gene', type=int, required=False, default=50)
parser.add_argument("--custom_genes", help="Additional genes to include", type=str, nargs="*")

args = parser.parse_args()

datasets = [ad.read_h5ad(f) for f in args.input]

adata = ad.concat(datasets)
adata_outer = ad.concat(datasets, join='outer')

additional_genes = [gene for gene in args.custom_genes if gene not in adata.var_names and gene in adata_outer.var_names]

# Add custom genes from outer join to the intersection
if additional_genes:
    adata_additional = adata_outer[adata.obs_names, additional_genes]
    adata_concatenated = ad.concat([adata, adata_additional], join="outer", axis=1)
    adata_concatenated.obs, adata_concatenated.obsm = adata.obs, adata.obsm
    adata = adata_concatenated

# Convert to CSR matrix
adata.X = csr_matrix(adata.X)
adata_outer.X = csr_matrix(adata_outer.X)

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

for column in adata.obs.columns:
    if column == "transfer":
        continue
    if not adata.obs[column].dtype.name == "category" and not adata.obs[column].dtype.name == "object":
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