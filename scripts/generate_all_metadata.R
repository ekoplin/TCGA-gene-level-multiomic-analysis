#!/usr/bin/env Rscript
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
script_files <- c(
  "generate_ensembl_gene_map.R",
  "generate_mirdb_targets.R",
  "generate_promoter_probes.R"
)

for (script_name in script_files) {
  source(file.path(script_dir, script_name))
}

cat("All metadata tables generated successfully!\n")
