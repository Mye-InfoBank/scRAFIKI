#!/usr/bin/env python

import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Integrate data using Harmony')
parser.add_argument('--input', help='Input file')
parser.add_argument('--output', help='Output file')
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

sc.pp.pca(adata)
sc.external.pp.harmony_integrate(adata, key='batch')

adata.write_h5ad(args.output)