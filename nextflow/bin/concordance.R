#!/usr/bin/env Rscript
# usage: concordance.R <primary_results_all.tsv> <secondary_results_all.tsv> <padj> <lfc>
a <- commandArgs(trailingOnly = TRUE); pc <- as.numeric(a[3]); lc <- as.numeric(a[4])
p <- read.delim(a[1], stringsAsFactors = FALSE); s <- read.delim(a[2], stringsAsFactors = FALSE)
sg <- function(d) d$gene_id[!is.na(d$padj) & d$padj < pc & abs(d$log2FoldChange) >= lc]
sp <- sg(p); ss <- sg(s); com <- intersect(p$gene_id, s$gene_id)
m <- merge(p[, c("gene_id", "log2FoldChange")], s[, c("gene_id", "log2FoldChange")], by = "gene_id", suffixes = c("_fc", "_salmon"))
m <- m[is.finite(m$log2FoldChange_fc) & is.finite(m$log2FoldChange_salmon), ]
out <- data.frame(metric = c("genes_with_padj_featurecounts", "genes_with_padj_salmon", "genes_in_both", "sig_featurecounts", "sig_salmon", "sig_overlap", "sig_union", "jaccard_sig",
                             "pearson_log2FC_all_common", "spearman_log2FC_all_common", "pearson_log2FC_sig_in_either"),
  value = c(sum(!is.na(p$padj)), sum(!is.na(s$padj)), length(com), length(sp), length(ss), length(intersect(sp, ss)), length(union(sp, ss)),
            round(length(intersect(sp, ss)) / max(1, length(union(sp, ss))), 4), round(cor(m[[2]], m[[3]]), 4), round(cor(m[[2]], m[[3]], method = "spearman"), 4),
            { e <- m[m$gene_id %in% union(sp, ss), ]; if (nrow(e) > 2) round(cor(e[[2]], e[[3]]), 4) else NA }))
write.table(out, "concordance_featurecounts_vs_salmon.tsv", sep = "\t", quote = FALSE, row.names = FALSE); print(out)
