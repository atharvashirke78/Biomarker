#!/usr/bin/env Rscript
# CHANGED vs copy of /workspace/biomarker/pipeline/bin/tximport.R: countsFromAbundance = "no", writes TPM.
# usage: tximport.R <tx2gene.tsv (no header: transcript_id, gene_id)> <dir holding <label>_quant/> <sample_map.tsv>   (sample IDs in the output = SRR)
suppressPackageStartupMessages({ library(tximport); library(readr) })
a <- commandArgs(trailingOnly = TRUE)
tx2gene <- read_tsv(a[1], col_names = c("transcript_id", "gene_id"), show_col_types = FALSE)
tx2gene <- tx2gene[!duplicated(tx2gene$transcript_id), ]
dirs <- sort(list.dirs(a[2], recursive = FALSE, full.names = FALSE)); dirs <- dirs[grepl("_quant$", dirs)]
if (length(dirs) == 0) stop("No *_quant directories found in ", a[2])
files <- file.path(a[2], dirs, "quant.sf"); smap <- read.delim(a[3], stringsAsFactors = FALSE)
lab <- sub("_quant$", "", dirs); srr <- smap$srr[match(lab, smap$label)]
if (anyNA(srr)) stop("No SRR in sample_map for label(s): ", paste(lab[is.na(srr)], collapse = ", "))
names(files) <- srr
stopifnot(all(file.exists(files)))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, countsFromAbundance = "no")
w <- function(m, f) write_tsv(data.frame(gene_id = rownames(m), m, check.names = FALSE), f)
w(txi$counts, "counts_salmon_tximport.tsv"); w(txi$abundance, "salmon_tpm.tsv"); saveRDS(txi, "salmon_tximport.rds")
cat(sprintf("tximport: %d genes x %d samples\n", nrow(txi$counts), ncol(txi$counts)))
