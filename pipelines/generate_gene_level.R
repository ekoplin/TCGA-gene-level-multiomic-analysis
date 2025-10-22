#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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

clean_header_row <- function(df, expected) {
  if (nrow(df) == 0) {
    return(df)
  }
  header_candidates <- tolower(as.character(df[1, seq_len(min(ncol(df), length(expected)))]))
  if (all(header_candidates == tolower(expected))) {
    df[-1, , drop = FALSE]
  } else {
    df
  }
}

make_mir_gcentric <- function(mir_file, expr_reference, mirna_targets) {
  mir_df <- utils::read.csv(mir_file, check.names = FALSE, stringsAsFactors = FALSE)
  mir_names <- mir_df[[1]]
  mir_mat <- as.matrix(mir_df[-1])
  rownames(mir_mat) <- mir_names
  colnames(mir_mat) <- colnames(mir_df)[-1]
  storage.mode(mir_mat) <- "numeric"

  samples <- intersect(colnames(mir_mat), colnames(expr_reference))
  if (length(samples) == 0) {
    stop("No shared samples between miRNA and expression matrices.", call. = FALSE)
  }
  mir_mat <- mir_mat[, samples, drop = FALSE]
  expr_ref <- expr_reference[, samples, drop = FALSE]

  targets_df <- utils::read.delim(mirna_targets, header = FALSE, stringsAsFactors = FALSE)
  targets_df <- clean_header_row(targets_df, c("mirna", "gene"))
  if (ncol(targets_df) < 2) {
    stop("miRNA target mapping must have at least two columns (miRNA and gene).", call. = FALSE)
  }
  gene_targets <- split(targets_df[[1]], targets_df[[2]])

  gene_names <- rownames(expr_ref)
  result <- matrix(0, nrow = length(gene_names), ncol = length(samples),
                   dimnames = list(gene_names, samples))

  for (gene in gene_names) {
    mi_list <- gene_targets[[gene]]
    if (length(mi_list) == 0) {
      next
    }
    mi_list <- intersect(mi_list, rownames(mir_mat))
    if (length(mi_list) == 0) {
      next
    }
    vals <- mir_mat[mi_list, , drop = FALSE]
    averaged <- colMeans(vals, na.rm = TRUE)
    averaged[is.nan(averaged)] <- 0
    result[gene, ] <- averaged
  }

  result
}

make_methylation_gcentric <- function(meth_file, expr_reference, transcript_map_file, promoter_file) {
  meth_df <- utils::read.csv(meth_file, check.names = FALSE, stringsAsFactors = FALSE)
  probe_ids <- meth_df[[1]]
  meth_mat <- as.matrix(meth_df[-1])
  rownames(meth_mat) <- probe_ids
  colnames(meth_mat) <- colnames(meth_df)[-1]
  storage.mode(meth_mat) <- "numeric"

  samples <- intersect(colnames(meth_mat), colnames(expr_reference))
  if (length(samples) == 0) {
    stop("No shared samples between methylation and expression matrices.", call. = FALSE)
  }
  meth_mat <- meth_mat[, samples, drop = FALSE]
  expr_ref <- expr_reference[, samples, drop = FALSE]

  transcript_df <- utils::read.delim(transcript_map_file, header = FALSE, stringsAsFactors = FALSE)
  transcript_df <- clean_header_row(transcript_df, c("gene", "transcript"))
  if (ncol(transcript_df) < 2) {
    stop("Transcript map must have at least two columns (gene and transcript).", call. = FALSE)
  }
  transcript_to_gene <- stats::setNames(transcript_df[[1]], transcript_df[[2]])

  promoter_df <- utils::read.delim(promoter_file, header = FALSE, stringsAsFactors = FALSE)
  promoter_df <- clean_header_row(promoter_df, c("probe", "chromosome", "start", "end", "transcript"))
  if (ncol(promoter_df) < 5) {
    stop("Promoter annotation must have at least five columns (probe ... transcript).", call. = FALSE)
  }
  gene_ids <- unname(transcript_to_gene[promoter_df[[5]]])
  valid <- !is.na(gene_ids)
  if (!any(valid)) {
    stop("No promoter probes map to genes using the provided transcript map.", call. = FALSE)
  }
  probe_by_gene <- split(promoter_df[[1]][valid], gene_ids[valid])

  gene_names <- rownames(expr_ref)
  result <- matrix(0, nrow = length(gene_names), ncol = length(samples),
                   dimnames = list(gene_names, samples))

  for (gene in gene_names) {
    probes <- probe_by_gene[[gene]]
    if (length(probes) == 0) {
      next
    }
    probes <- intersect(probes, rownames(meth_mat))
    if (length(probes) == 0) {
      next
    }
    vals <- meth_mat[probes, , drop = FALSE]
    averaged <- colMeans(vals, na.rm = TRUE)
    averaged[is.nan(averaged)] <- 0
    result[gene, ] <- averaged
  }

  result
}

expr_mat <- read_matrix(rna_file)

message("Generating methylation gene-level matrix...")
meth_mat <- make_methylation_gcentric(meth_file, expr_mat, transcript_map, promoter_file)

message("Generating miRNA gene-level matrix...")
mir_mat <- make_mir_gcentric(mirna_file, expr_mat, mirna_mapping)

message("Preparing gene expression matrix...")
expr_sub <- expr_mat

common_samples <- Reduce(intersect, list(colnames(expr_sub), colnames(meth_mat), colnames(mir_mat)))
if (length(common_samples) == 0) {
  stop("No shared samples across omics layers after harmonisation.", call. = FALSE)
}
expr_sub <- expr_sub[, common_samples, drop = FALSE]
meth_mat <- meth_mat[, common_samples, drop = FALSE]
mir_mat <- mir_mat[, common_samples, drop = FALSE]

expr_norm <- normalize_matrix(expr_sub)
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
