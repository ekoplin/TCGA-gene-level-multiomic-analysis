#!/usr/bin/env Rscript
library(data.table)

input_path <- "metadata/miRDB_v6.0_prediction_result_hsa.txt.gz"

# Read all columns (robust to different headers)
mir <- fread(input_path)
cols <- tolower(names(mir))

# Normalize header names
if ("mirna" %in% cols) setnames(mir, "miRNA", "miRNA_Name")
if ("target_gene" %in% cols) setnames(mir, "Target_Gene", "Target_Gene")
if ("target gene" %in% cols) setnames(mir, "Target Gene", "Target_Gene")

df <- unique(mir[, .(miRNA_Name, Target_Gene)])
write.table(df, "metadata/mirna_target_mirdb.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("✅ metadata/mirna_target_mirdb.txt generated successfully")
