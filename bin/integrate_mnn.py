#!/usr/bin/env python3

import mnnpy
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Integrate data using MNN')
parser.add_argument('--input', help='Input file', type=str)
parser.add_argument('--output', help='Output file', type=str)
parser.add_argument('--cpus', help='Number of CPUs to use', default=8, type=int)

args = parser.parse_args()

adata = sc.read_h5ad(args.input)

batches = adata.obs["batch"].unique()

# split by batch for mnnpy
adata_list = [adata[adata.obs["batch"] == batch] for batch in batches]

# integration with mnn
corrected = mnnpy.mnn_correct(*adata_list, batch_categories=batches, n_jobs=args.cpus)

corrected_adata = corrected[0]
corrected_adata.write_h5ad(args.output)