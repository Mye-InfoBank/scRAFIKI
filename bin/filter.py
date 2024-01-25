#!/opt/conda/bin/python

import scanpy as sc
import numpy as np
import argparse
import mygene

parser = argparse.ArgumentParser(description="Filter dataset")
parser.add_argument("--input", help="Input file", type=str)
parser.add_argument("--id", help="Dataset ID", type=str)
parser.add_argument("--output", help="Output file", type=str)

parser.add_argument("--min_genes", help="Minimum number of genes", type=int, required=False)
parser.add_argument("--max_genes", help="Maximum number of genes", type=int, required=False)
parser.add_argument("--min_counts", help="Minimum number of counts", type=int, required=False)
parser.add_argument("--max_counts", help="Maximum number of counts", type=int, required=False)
parser.add_argument("--max_pct_mito", help="Maximum percentage of mitochondrial counts", type=float, required=False)

args = parser.parse_args()

# Function borrowed from https://github.com/icbi-lab/luca/blob/5ffb0a4671e9c288b10e73de18d447ee176bef1d/lib/scanpy_helper_submodule/scanpy_helpers/util.py#L122C1-L135C21
def aggregate_duplicate_var(adata, aggr_fun=np.mean):
    retain_var = ~adata.var_names.duplicated(keep="first")
    duplicated_var = adata.var_names[adata.var_names.duplicated()].unique()
    if len(duplicated_var):
        for var in duplicated_var:
            mask = adata.var_names == var
            var_aggr = aggr_fun(adata.X[:, mask], axis=1)[:, np.newaxis]
            adata.X[:, mask] = np.repeat(var_aggr, np.sum(mask), axis=1)

        adata_dedup = adata[:, retain_var].copy()
        return adata_dedup
    else:
        return adata

adata = sc.read_h5ad(args.input)
adata.obs["dataset"] = args.id

if adata.__dict__["_raw"] and "_index" in adata.__dict__["_raw"].__dict__["_var"]:
    adata.__dict__["_raw"].__dict__["_var"] = (
        adata.__dict__["_raw"].__dict__["_var"].rename(columns={"_index": "features"})
    )

# Convert everything to gene symbols
mg = mygene.MyGeneInfo()
df_genes = mg.querymany(adata.var.index, scopes=["symbol", "entrezgene", "ensemblgene"], fields="symbol", species="human", as_dataframe=True) 
mapping = df_genes["symbol"].dropna().to_dict()

adata.var_names = adata.var.index.map(lambda x: mapping.get(x, x))

# Convert varnames to upper case
adata.var_names = adata.var_names.str.upper()
adata.var_names = adata.var_names.str.replace("_", "-")

# Calculate mean of same-named genes
adata = aggregate_duplicate_var(adata)

if "mito" not in adata.var.columns:
    adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
)

# very basic gene filtering - genes with 0 cells cause some downstream processes to fail.
print("Filtering genes")
print(f"    Before: {adata.shape[1]}")
sc.pp.filter_genes(adata, min_counts=3)
print(f"    After: {adata.shape[1]}")

# Apply thresholds
if args.min_counts:
    print("Filter by min_counts")
    print(f"    Before: {adata.shape[0]}")
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    print(f"    After: {adata.shape[0]}")

if args.max_counts:
    print("Filter by max_counts")
    print(f"    Before: {adata.shape[0]}")
    sc.pp.filter_cells(adata, max_counts=args.max_counts)
    print(f"    After: {adata.shape[0]}")

if args.min_genes:
    print("Filter by min_genes")
    print(f"    Before: {adata.shape[0]}")
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    print(f"    After: {adata.shape[0]}")

if args.max_genes:
    print("Filter by max_genes")
    print(f"    Before: {adata.shape[0]}")
    sc.pp.filter_cells(adata, max_genes=args.max_genes)
    print(f"    After: {adata.shape[0]}")

if args.max_pct_mito:
    print("Filter by max_pct_mito")
    print(f"    Before: {adata.shape[0]}")
    adata = adata[adata.obs["pct_counts_mito"] < args.max_pct_mito].copy()
    print(f"    After: {adata.shape[0]}")

# Print all non-unique varnames
non_unique_varnames = adata.var_names[adata.var_names.duplicated(keep=False)]

if len(non_unique_varnames) > 0:
    raise ValueError(f"Non-unique varnames found: {non_unique_varnames}")

adata.write_h5ad(args.output)
