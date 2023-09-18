#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("At least two arguments must be supplied (annData file, batch)",
    call. = FALSE)
}

library(anndata)
library(celda)

batch <- args[2]
ad <- read_h5ad(args[1])

ad <- ad[ad$obs$batch == batch, ]$copy()

results <- decontX(ad$T$X)

ad$X <- t(results$decontXcounts)

write_h5ad(ad, paste0(batch, "_ambient.h5ad"))
