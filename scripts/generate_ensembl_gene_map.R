#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
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

gtf_path <- file.path(metadata_dir, "Homo_sapiens.GRCh38.96.gtf.gz")
if (!file.exists(gtf_path)) {
  stop("Input GTF not found: ", gtf_path)
}

gtf <- import(gtf_path)
transcripts <- gtf[gtf$type == "transcript"]
if (length(transcripts) == 0) {
  stop("No transcript features found in GTF.")
}

transcript_metadata <- mcols(transcripts)
required_cols <- c("transcript_id", "gene_id", "gene_name")
missing_cols <- setdiff(required_cols, colnames(transcript_metadata))
if (length(missing_cols) > 0) {
  stop("Missing required columns in GTF metadata: ", paste(missing_cols, collapse = ", "))
}

transcript_df <- as.data.frame(transcript_metadata[required_cols])
colnames(transcript_df) <- c("transcript_id", "gene_id", "gene_name")
transcript_df <- unique(transcript_df)

output_path <- file.path(metadata_dir, "ensembl_gene_transcript_map.txt")
write.table(
  transcript_df,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Ensembl gene-transcript map generated at ", output_path)
