#!/opt/conda/bin/python

import scanpy as sc
import numpy as np
import argparse
from scipy.sparse import csr_matrix

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

parser = argparse.ArgumentParser(description="Filter dataset")
parser.add_argument("--input", help="Input file", type=str)
parser.add_argument("--id", help="Dataset ID", type=str)
parser.add_argument("--output", help="Output file", type=str)
parser.add_argument("--problems", help="Problems file", type=str)
parser.add_argument("--no-symbols", help="Convert varnames to gene symbols", action="store_true")
parser.add_argument("--transfer", help="Apply transfer leanring on dataset", action="store_true")
parser.add_argument("--custom_metadata", help="Additional metadata columns to include", type=str, nargs="*")

parser.add_argument("--min_genes", help="Minimum number of genes", type=int, required=False)
parser.add_argument("--max_genes", help="Maximum number of genes", type=int, required=False)
parser.add_argument("--min_counts", help="Minimum number of counts", type=int, required=False)
parser.add_argument("--max_counts", help="Maximum number of counts", type=int, required=False)
parser.add_argument("--max_pct_mito", help="Maximum percentage of mitochondrial counts", type=float, required=False)

args = parser.parse_args()

problems = []

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
    
columns_additional = {column: False for column in args.custom_metadata}
columns_considered = {**columns_required, **columns_additional}

print("Reading input")
adata = sc.read_h5ad(args.input)
adata.obs["dataset"] = args.id
adata.obs["transfer"] = args.transfer

if adata.__dict__["_raw"] and "_index" in adata.__dict__["_raw"].__dict__["_var"]:
    adata.__dict__["_raw"].__dict__["_var"] = (
        adata.__dict__["_raw"].__dict__["_var"].rename(columns={"_index": "features"})
    )

# Make sure adata.X contains raw counts
max_diff = np.abs(adata.X - np.round(adata.X)).max()
if max_diff > 1e-3:
    problems.append(f"adata.X does not contain raw counts. Max diff: {max_diff}")

# Make sure dataset is not empty
if adata.shape[0] == 0 or adata.shape[1] == 0:
    problems.append(f"Dataset is empty.")

# Make sure all required columns are present
for column, required in columns_considered.items():
    if column not in adata.obs.columns:
        if required:
            problems.append(f"Column {column} is required but not found.")
        else:
            adata.obs[column] = "Unknown"

if not args.no_symbols:
    # Make sure there are no ensembl gene ids or entrez gene ids in varnames
    if adata.var_names.str.match(r'ENS[A-Z]+[0-9]{11}').any():
        problems.append("varnames contain ensembl gene ids")
    if adata.var_names.str.match(r"d+").any():
        problems.append("varnames contain entrez gene ids (or just numbers)")

# Make sure there are no duplicate obs names
if adata.obs_names.duplicated().any():
    problems.append("Duplicate obs names found")

# Make sure there are no duplicate var names
if adata.var_names.duplicated().any():
    problems.append("Duplicate var names found")

# Terminate if there are any problems
if problems:
    with open(args.problems, "w") as f:
        f.write("\n".join(problems))
    print("Problems found. Terminating.")
    exit(0)

# Subset columns
adata.obs = adata.obs[columns_considered.keys()]

# Prepare obs names
adata.obs_names = args.id + "_" + adata.obs_names
adata.obs_names = adata.obs_names.str.replace("_", "-")

# Convert to CSR matrix
adata.X = csr_matrix(adata.X)

if args.no_symbols:
    print("Converting varnames to gene symbols")
    import mygene
    mg = mygene.MyGeneInfo()
    df_genes = mg.querymany(adata.var.index, scopes=["symbol", "entrezgene", "ensemblgene"], fields="symbol", species="human", as_dataframe=True) 
    mapping = df_genes["symbol"].dropna().to_dict()

    adata.var_names = adata.var.index.map(lambda x: mapping.get(x, x))

# Convert varnames to upper case
adata.var_names = adata.var_names.str.upper()
adata.var_names = adata.var_names.str.replace("_", "-")
adata.var_names = adata.var_names.str.replace(".", "-")

print("Aggregating duplicate varnames")
# Calculate mean of same-named genes
adata = aggregate_duplicate_var(adata)

print("Calculating QC metrics")
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

print("Writing output")
adata.write_h5ad(args.output)
