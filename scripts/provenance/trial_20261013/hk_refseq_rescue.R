#!/usr/bin/env Rscript
# Annotation lookup only (no new install): RefSeq NM accession -> Entrez -> Ensembl gene ID with org.Hs.eg.db, for HK-list genes whose symbol is NOT found in GENCODE v50.
# Run with env/prefix/r_pathway/bin/Rscript (the only env that already has org.Hs.eg.db); the rest of the analysis runs in r_dge.
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
hk <- read.table("results/13_reference_genes/source/HK_genes.txt", sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, strip.white = TRUE, col.names = c("symbol", "refseq"), fill = TRUE)
gs <- read.delim("ref/gene_id_to_symbol.tsv", stringsAsFactors = FALSE)
un <- hk[!(hk$symbol %in% gs$gene_name), ]
acc <- sub("\\..*$", "", un$refseq); ok <- nzchar(acc)
m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(acc[ok]), keytype = "REFSEQ", columns = c("ENTREZID", "ENSEMBL", "SYMBOL")))
m <- m[!is.na(m$ENSEMBL), ]
r <- merge(data.frame(symbol = un$symbol, refseq = un$refseq, REFSEQ = acc, stringsAsFactors = FALSE), m, by = "REFSEQ")
write.table(r[, c("symbol", "refseq", "ENTREZID", "SYMBOL", "ENSEMBL")], "results/13_reference_genes/hk_refseq_rescue_map.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("HK symbols not in GENCODE v50:", nrow(un), "| with a RefSeq accession:", sum(ok), "| rescued to >=1 Ensembl ID:", length(unique(r$symbol)), "\n")
writeLines(c(paste("org.Hs.eg.db", packageVersion("org.Hs.eg.db")), paste("R", getRversion())), "results/13_reference_genes/tool_versions_rescue.txt")
