#!/usr/bin/env Rscript

# This file is part of the scib-pipeline: https://github.com/theislab/scib-pipeline/blob/main/scripts/integration/runMethods.R

library('optparse')
library(rlang)
require(Seurat)

methods <- c('seurat', 'seuratrpca', 'conos', 'harmony', 'liger', 'fastmnn')

option_list <- list(make_option(c("-m", "--method"), type="character", default=NA, help="integration method to use", choices=methods),
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
		    make_option(c("-b", "--batch"), type="character", default=NA, help="batch variable"),
		    make_option(c("-v", "--hvg"), type="character", default=NA, help="hvg list for seurat"))



opt = parse_args(OptionParser(option_list=option_list))

source('integration.R')
sobj = loadSeuratObject(opt$i)

if(opt$method=='seurat'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	else {
		hvg <- rownames(sobj@assays$RNA)
	}

	out = runSeurat(sobj, opt$b, hvg)
}

if(opt$method=='seuratrpca'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	else {
		hvg <- rownames(sobj@assays$RNA)
	}

	out = runSeuratRPCA(sobj, opt$b, hvg)
}


if(opt$method=='conos'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
		sobj <- subset(sobj, features=hvg)
	}
	out = runConos(sobj, opt$b)
}

if(opt$method=='harmony'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
		sobj <- subset(sobj, features=hvg)
	}

	out=runHarm(sobj, opt$b)
}

if(opt$method=='liger'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	else {
		hvg <- rownames(sobj@assays$RNA)
	}

	out = runLiger(sobj, opt$b, hvg)
}

if(opt$method=='fastmnn'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
		sobj <- subset(sobj, features=hvg)
	}

	out=runFastMNN(sobj, opt$b)
}

saveSeuratObject(out, opt$o)
