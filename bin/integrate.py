#!/usr/bin/env python3

import anndata as ad
import scib
import argparse

# Disable warnings
import warnings
warnings.filterwarnings("ignore")

methods = {
    "bbknn": scib.ig.bbknn,
    "combat": scib.ig.combat,
    "desc": scib.ig.desc,
    "harmony": scib.ig.harmony,
    "mnn": scib.ig.mnn,
    #"saucie": scib.ig.saucie, # Deactivated because dependencies are not available
    "scanorama": scib.ig.scanorama,
    "scanvi": scib.ig.scanvi,
    "scgen": scib.ig.scgen,
    "scvi": scib.ig.scvi,
    #"trvae": scib.ig.trvae, # Deactivated because it requires keras==2.2.4 which is incompatible with other dependencies
    "trvaep": scib.ig.trvaep,
}

parser = argparse.ArgumentParser(description='Integrate data')
parser.add_argument('--input', help='Input file', type=str)
parser.add_argument('--output', help='Output file', type=str)
parser.add_argument('--method', help='Integration method', type=str, choices=methods.keys())
parser.add_argument('--scvi_model', help='scvi model', type=str)

args = parser.parse_args()

adata = ad.read_h5ad(args.input)

method = methods[args.method]

if args.method == "scanvi":
    if args.scvi_model is None:
        adata = method(adata, "batch", "celltype", save_model=True)
    else:
        adata = method(adata, "batch", "celltype", scvi_model_path=args.scvi_model, save_model=True)
elif args.method == "scgen":
        adata = method(adata, "batch", "celltype")
elif args.method == "scvi":
    adata = method(adata, "batch", save_model=True)
else:
    adata = method(adata, "batch")

adata.write_h5ad(args.output)