#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
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

input_path <- file.path(metadata_dir, "miRDB_v6.0_prediction_result.txt.gz")
if (!file.exists(input_path)) {
  stop("miRDB prediction file not found: ", input_path)
}

dt <- fread(input_path, select = c("miRNA_Name", "Target_Gene"))
if (nrow(dt) == 0) {
  stop("No records found in miRDB file.")
}

dt_unique <- unique(dt)

output_path <- file.path(metadata_dir, "mirna_target_mirdb.txt")
fwrite(dt_unique, file = output_path, sep = "\t", quote = FALSE)

message("miRDB target table generated at ", output_path)
