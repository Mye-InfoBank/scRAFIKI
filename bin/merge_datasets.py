#!/opt/conda/bin/python

import argparse
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix

parser = argparse.ArgumentParser(description="Merge datasets")
parser.add_argument("--input", help="Input file", type=str, nargs="+")
parser.add_argument("--base", help="Base dataset to use as reference", type=str, required=False)
parser.add_argument("--output_intersection", help="Output file containing all cells but gene intersection", type=str, required=True)
parser.add_argument("--output_union", help="Output file, outer join of cells and genes", type=str, required=True)
parser.add_argument("--output_transfer", help="Output file, cells to project onto base", type=str, required=False)
parser.add_argument("--min_cells", help='Minimum number of cells to keep a gene', type=int, required=False, default=50)
parser.add_argument("--custom_genes", help="Additional genes to include", type=str, nargs="*")

args = parser.parse_args()

datasets = [ad.read_h5ad(f) for f in args.input]

if args.base:
    if not args.output_transfer:
        raise ValueError("Transfer file required when using base dataset")

    adata_base = ad.read_h5ad(args.base)
    datasets = [adata_base] + datasets

adata_intersection = ad.concat(datasets)
adata_union = ad.concat(datasets, join='outer')

additional_genes = [gene for gene in args.custom_genes if gene not in adata_intersection.var_names and gene in adata_union.var_names]

# Add custom genes from outer join to the intersection
if additional_genes:
    adata_additional = adata_union[adata_intersection.obs_names, additional_genes]
    adata_concatenated = ad.concat([adata_intersection, adata_additional], join="outer", axis=1)
    adata_concatenated.obs, adata_concatenated.obsm = adata_intersection.obs, adata_intersection.obsm
    adata_intersection = adata_concatenated

# Convert to CSR matrix
adata_intersection.X = csr_matrix(adata_intersection.X)
adata_union.X = csr_matrix(adata_union.X)

# Filter cells with no counts
cell_mask, _ = sc.pp.filter_cells(adata_intersection, min_genes=1, inplace=False)
adata_intersection = adata_intersection[cell_mask, :]
adata_union = adata_union[cell_mask, :]

# Filter genes with too few occurrences in outer join
sc.pp.filter_genes(adata_union, min_cells=args.min_cells)

adata_intersection.obs["batch"] = adata_intersection.obs["dataset"].astype(str) + "_" + adata_intersection.obs["batch"].astype(str)
adata_intersection.obs["patient"] = adata_intersection.obs["dataset"].astype(str) + "_" + adata_intersection.obs["patient"].astype(str)

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

for column in adata_intersection.obs.columns:
    if not adata_intersection.obs[column].dtype.name == "category" and not adata_intersection.obs[column].dtype.name == "object":
        continue
    # Convert first to string and then to category
    adata_intersection.obs[column] = adata_intersection.obs[column].astype(str).fillna("Unknown").apply(to_Florent_case).astype("category")

adata_union.obs = adata_intersection.obs

adata_intersection.layers["counts"] = adata_intersection.X
adata_union.layers["counts"] = adata_union.X

adata_intersection.write_h5ad(args.output_intersection)
adata_union.write_h5ad(args.output_union)

if args.base:
    adata_transfer = adata_intersection[~adata_intersection.obs.index.isin(adata_base.obs.index)]
    adata_transfer.write_h5ad(args.output_transfer)