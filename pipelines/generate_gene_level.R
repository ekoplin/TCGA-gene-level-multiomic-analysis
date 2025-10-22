#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(reticulate)
  library(yaml)
  library(preprocessCore)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path, call. = FALSE)
}

config <- yaml::read_yaml(config_path)
project <- config$project %||% "TCGA-BRCA"
data_dir <- config$data_dir %||% "data"
output_dir <- config$output_dir %||% "output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- config$inputs %||% list()
rna_file <- inputs$rna %||% file.path(data_dir, paste0(project, "_rna_expression.csv"))
meth_file <- inputs$methylation %||% file.path(data_dir, paste0(project, "_dna_methylation.csv"))
mirna_file <- inputs$mirna %||% file.path(data_dir, paste0(project, "_mirna_expression.csv"))

metadata <- config$metadata %||% list()
mirna_mapping <- metadata$mirna_targets
promoter_file <- metadata$methylation_promoters
transcript_map <- metadata$gene_transcript_map

if (any(vapply(list(mirna_mapping, promoter_file, transcript_map), is.null, logical(1)))) {
  stop("Configuration must provide metadata paths: mirna_targets, methylation_promoters, and gene_transcript_map.", call. = FALSE)
}

data_inputs <- c(rna_file, meth_file, mirna_file)
missing_data <- data_inputs[!file.exists(data_inputs)]
if (length(missing_data) > 0) {
  stop("Missing required omics matrices:\n", paste(" -", missing_data, collapse = "\n"), call. = FALSE)
}

metadata_inputs <- c(mirna_mapping, promoter_file, transcript_map)
missing_meta <- metadata_inputs[!file.exists(metadata_inputs)]
if (length(missing_meta) > 0) {
  stop("Missing required metadata files:\n", paste(" -", missing_meta, collapse = "\n"), call. = FALSE)
}

message("Loading python MONTI utilities...")
monti <- reticulate::import_from_path("monti", path = "src")

read_matrix <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- df[[1]]
  mat <- as.matrix(df[-1])
  rownames(mat) <- ids
  colnames(mat) <- colnames(df)[-1]
  storage.mode(mat) <- "numeric"
  mat
}

normalize_matrix <- function(mat) {
  storage.mode(mat) <- "numeric"
  mat <- log2(mat + 1)
  qn <- preprocessCore::normalize.quantiles(mat)
  rownames(qn) <- rownames(mat)
  colnames(qn) <- colnames(mat)
  scaled <- t(apply(qn, 1, function(x) {
    rng <- max(x) - min(x)
    if (is.na(rng) || rng == 0) {
      rep(0, length(x))
    } else {
      (x - min(x)) / rng
    }
  }))
  rownames(scaled) <- rownames(mat)
  colnames(scaled) <- colnames(mat)
  scaled
}

message("Generating methylation gene-level matrix...")
meth_gene <- monti$make_methylation_gcentric(meth_file, rna_file, transcript_map, promoter_file)
meth_mat <- as.matrix(meth_gene)
storage.mode(meth_mat) <- "numeric"

message("Generating miRNA gene-level matrix...")
mir_gene <- monti$make_mir_gcentric(mirna_file, rna_file, mirna_mapping)
mir_mat <- as.matrix(mir_gene)
storage.mode(mir_mat) <- "numeric"

message("Preparing gene expression matrix...")
expr_mat <- read_matrix(rna_file)

common_samples <- Reduce(intersect, list(colnames(expr_mat), colnames(meth_mat), colnames(mir_mat)))
if (length(common_samples) == 0) {
  stop("No shared samples across omics layers after harmonisation.", call. = FALSE)
}
expr_mat <- expr_mat[, common_samples, drop = FALSE]
meth_mat <- meth_mat[, common_samples, drop = FALSE]
mir_mat <- mir_mat[, common_samples, drop = FALSE]

expr_norm <- normalize_matrix(expr_mat)
meth_norm <- normalize_matrix(meth_mat)
mir_norm <- normalize_matrix(mir_mat)

common_genes <- Reduce(intersect, list(rownames(expr_norm), rownames(meth_norm), rownames(mir_norm)))
if (length(common_genes) == 0) {
  stop("No shared genes across omics layers after harmonisation.", call. = FALSE)
}
common_genes <- sort(common_genes)
expr_norm <- expr_norm[common_genes, , drop = FALSE]
meth_norm <- meth_norm[common_genes, , drop = FALSE]
mir_norm <- mir_norm[common_genes, , drop = FALSE]

write_matrix <- function(mat, path) {
  df <- data.frame(Gene = rownames(mat), mat, check.names = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

expr_out <- file.path(output_dir, "gene_expression_matrix.tsv")
meth_out <- file.path(output_dir, "gene_methylation_matrix.tsv")
mir_out <- file.path(output_dir, "gene_miRNA_matrix.tsv")

write_matrix(expr_norm, expr_out)
write_matrix(meth_norm, meth_out)
write_matrix(mir_norm, mir_out)

message("Saved gene-level matrices to ", output_dir)
