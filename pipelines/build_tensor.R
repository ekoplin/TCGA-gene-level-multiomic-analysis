#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path, call. = FALSE)
}

config <- yaml::read_yaml(config_path)
output_dir <- config$output_dir %||% "output"
matrices <- config$matrices %||% list()

expr_path <- matrices$expression %||% file.path(output_dir, "gene_expression_matrix.tsv")
meth_path <- matrices$methylation %||% file.path(output_dir, "gene_methylation_matrix.tsv")
mir_path <- matrices$mirna %||% file.path(output_dir, "gene_miRNA_matrix.tsv")

load_gene_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Matrix not found: ", path, call. = FALSE)
  }
  df <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  genes <- df[[1]]
  mat <- as.matrix(df[-1])
  rownames(mat) <- genes
  colnames(mat) <- colnames(df)[-1]
  storage.mode(mat) <- "numeric"
  mat
}

expr <- load_gene_matrix(expr_path)
meth <- load_gene_matrix(meth_path)
mir <- load_gene_matrix(mir_path)

common_genes <- Reduce(intersect, list(rownames(expr), rownames(meth), rownames(mir)))
if (length(common_genes) == 0) {
  stop("No shared genes across matrices. Ensure gene-level preprocessing was successful.", call. = FALSE)
}
common_genes <- sort(common_genes)

common_samples <- Reduce(intersect, list(colnames(expr), colnames(meth), colnames(mir)))
if (length(common_samples) == 0) {
  stop("No shared samples across matrices. Ensure gene-level preprocessing was successful.", call. = FALSE)
}
common_samples <- sort(common_samples)

expr <- expr[common_genes, common_samples, drop = FALSE]
meth <- meth[common_genes, common_samples, drop = FALSE]
mir <- mir[common_genes, common_samples, drop = FALSE]

omics <- c("expression", "methylation", "mirna")
tensor <- array(0, dim = c(length(common_genes), length(common_samples), length(omics)),
                dimnames = list(common_genes, common_samples, omics))

tensor[, , "expression"] <- expr
tensor[, , "methylation"] <- meth
tensor[, , "mirna"] <- mir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(tensor, file = file.path(output_dir, "tensor.rds"))
write.table(common_genes, file = file.path(output_dir, "genelist.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

message("✅ Tensor successfully built and saved as output/tensor.rds")
