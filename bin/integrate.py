#!/usr/bin/env python3

import anndata as ad
import scib
import argparse

methods = {
    "bbknn": scib.ig.bbknn,
    "combat": scib.ig.combat,
    "desc": scib.ig.desc,
    "harmony": scib.ig.harmony,
    "mnn": scib.ig.mnn,
    #"saucie": scib.ig.saucie, # Deactivated because dependencies are not available
    "scanorama": scib.ig.scanorama,
    "scanvi": scib.ig.scanvi,
    #"scgen": scib.ig.scgen, # Deactivated because of problems with deprecated dependencies
    "scvi": scib.ig.scvi,
    #"trvae": scib.ig.trvae, # Deactivated because it requires keras==2.2.4 which is incompatible with other dependencies
    "trvaep": scib.ig.trvaep,
}

parser = argparse.ArgumentParser(description='Integrate data')
parser.add_argument('--input', help='Input file', type=str)
parser.add_argument('--output', help='Output file', type=str)
parser.add_argument('--method', help='Integration method', type=str, choices=methods.keys())

args = parser.parse_args()

adata = ad.read_h5ad(args.input)

method = methods[args.method]

if args.method in ["scgen", "scanvi"]:
    adata = method(adata, "batch", "cell_type")
else:
    adata = method(adata, "batch")

adata.write_h5ad(args.output)