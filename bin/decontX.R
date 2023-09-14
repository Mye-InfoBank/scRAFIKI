args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(anndata)
library(MatrixExtra)
library(celda)

ad <- read_h5ad(args[1])
counts <- ad$T$X
counts <- as.csc.matrix(counts)

results <- decontX(counts)

decontXcounts <- results$decontXcounts

ad$layers[["before_ambient_correction"]] <- ad$X
ad$X <- t(decontXcounts)

write_h5ad(ad, "after_ambient_correction.h5ad")

