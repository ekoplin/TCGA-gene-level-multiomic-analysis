#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(yaml)
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
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

trim_barcodes <- function(barcodes) {
  substr(barcodes, 1, 16)
}

write_matrix <- function(mat, path, id_label, strip_version = FALSE) {
  mat <- as.matrix(mat)
  if (strip_version) {
    rownames(mat) <- sub("\\..*$", "", rownames(mat))
  }
  if (anyDuplicated(rownames(mat))) {
    message("Collapsing duplicated identifiers in ", basename(path), " by averaging.")
    idx <- rownames(mat)
    sums <- rowsum(mat, idx)
    ones <- matrix(1, nrow = nrow(mat), ncol = ncol(mat))
    colnames(ones) <- colnames(mat)
    counts <- rowsum(ones, idx)
    mat <- sums / counts
  }
  df <- data.frame(id = rownames(mat), mat, check.names = FALSE)
  colnames(df)[1] <- id_label
  utils::write.csv(df, path, row.names = FALSE)
}

fetch_and_save <- function(query_params, prefix, id_label, strip_version = FALSE) {
  message("Querying ", prefix, " data for ", project, "...")
  query <- do.call(GDCquery, c(list(project = project), query_params))
  GDCdownload(query)
  se <- GDCprepare(query)
  mat <- SummarizedExperiment::assay(se)
  colnames(mat) <- make.unique(trim_barcodes(colnames(mat)))
  outfile <- file.path(data_dir, paste0(project, "_", prefix, ".csv"))
  write_matrix(mat, outfile, id_label, strip_version = strip_version)
  saveRDS(se, file.path(data_dir, paste0(project, "_", prefix, ".rds")))
  message("Saved ", prefix, " matrix to ", outfile)
  invisible(outfile)
}

# RNA-seq
# rna_params <- config$download$rna
# if (is.null(rna_params)) {
#   rna_params <- list(
#     data.category = "Transcriptome Profiling",
#     data.type = "Gene Expression Quantification",
#     workflow.type = "HTSeq - FPKM"
#   )
# }
rna_params <- config$download$rna
if (is.null(rna_params)) {
  rna_params <- list(
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
}
fetch_and_save(rna_params, "rna_expression", "gene_id", strip_version = TRUE)

# DNA methylation
meth_params <- config$download$methylation
if (is.null(meth_params)) {
  meth_params <- list(
    data.category = "DNA Methylation",
    platform = "Illumina Human Methylation 450",
    legacy = TRUE
  )
}
fetch_and_save(meth_params, "dna_methylation", "probe_id")

# miRNA
mirna_params <- config$download$mirna
if (is.null(mirna_params)) {
  mirna_params <- list(
    data.category = "Transcriptome Profiling",
    data.type = "miRNA Expression Quantification",
    workflow.type = "BCGSC miRNA Profiling"
  )
}
fetch_and_save(mirna_params, "mirna_expression", "mirna_id")
