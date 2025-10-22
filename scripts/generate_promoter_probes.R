#!/usr/bin/env Rscript
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path)))
  }
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()
metadata_dir <- file.path(script_dir, "..", "metadata")
if (!dir.exists(metadata_dir)) {
  stop("metadata directory not found: ", metadata_dir)
}

data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_df <- as.data.frame(annotation)

required_cols <- c("Name", "UCSC_RefGene_Name", "Distance_to_TSS")
missing_cols <- setdiff(required_cols, colnames(annotation_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in annotation: ", paste(missing_cols, collapse = ", "))
}

annotation_df <- annotation_df[!is.na(annotation_df$Distance_to_TSS), ]
promoter_df <- annotation_df[abs(annotation_df$Distance_to_TSS) <= 2000, ]

result_df <- unique(data.frame(
  probe_ID = promoter_df$Name,
  gene_symbol = promoter_df$UCSC_RefGene_Name,
  stringsAsFactors = FALSE
))

output_path <- file.path(metadata_dir, "promoter_probes_illumina450.txt")
write.table(
  result_df,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Promoter probe table generated at ", output_path)
