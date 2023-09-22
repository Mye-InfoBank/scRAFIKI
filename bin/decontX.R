#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("At least one argument must be supplied (annData file)",
    call. = FALSE)
}

library(anndata)
library(celda)

ad <- read_h5ad(args[1])

results <- decontX(ad$T$X)

ad$X <- t(results$decontXcounts)

basename <- gsub(".h5ad", "", basename(args[1]))

write_h5ad(ad, paste0(basename, "_ambient.h5ad"))
