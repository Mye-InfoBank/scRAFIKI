#!/opt/conda/bin/python

import scanpy as sc
import argparse

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

adata = sc.read_h5ad(args.input)
adata.obs["dataset"] = args.id

if adata.__dict__["_raw"] and "_index" in adata.__dict__["_raw"].__dict__["_var"]:
    adata.__dict__["_raw"].__dict__["_var"] = (
        adata.__dict__["_raw"].__dict__["_var"].rename(columns={"_index": "features"})
    )

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

adata.write_h5ad(args.output)
