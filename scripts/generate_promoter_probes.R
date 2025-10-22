#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Required packages
for (pkg in c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
              "TxDb.Hsapiens.UCSC.hg19.knownGene", "GenomicRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Load probe annotation
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gr_probes <- GRanges(seqnames = anno$chr,
                     ranges = IRanges(start = anno$pos, width = 1),
                     strand = anno$strand)
mcols(gr_probes)$probe <- rownames(anno)

# Load gene TSS information
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tss <- promoters(genes(txdb), upstream = 0, downstream = 1)
tss_gr <- resize(tss, width = 1, fix = "start")

# Compute distance from probes to nearest TSS
overlap <- distanceToNearest(gr_probes, tss_gr)
distances <- rep(Inf, length(gr_probes))
distances[queryHits(overlap)] <- mcols(overlap)$distance

# Keep probes within ±2 kb of any TSS
promoter_idx <- which(distances <= 2000)
promoter_probes <- gr_probes[promoter_idx]

# Extract UCSC gene names if available
gene_names <- sapply(strsplit(as.character(anno$UCSC_RefGene_Name[promoter_idx]), ","), `[`, 1)

df <- data.frame(
  probe = promoter_probes$probe,
  gene = gene_names
)

output_path <- "metadata/promoter_probes_illumina450.txt"
write.table(df, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste0("✅ Promoter probe table generated at ", normalizePath(output_path)))
