#!/usr/bin/env Rscript
# Install core dependencies for the TCGA gene-level multi-omics pipeline

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("TCGAbiolinks", "minfi", "limma", "edgeR", "preprocessCore"), update = FALSE, ask = FALSE)

if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml")
}
