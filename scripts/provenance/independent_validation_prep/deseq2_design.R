#!/usr/bin/env Rscript
# NEW (adapted from nextflow/bin/deseq2.R, matrix route only; that file is NOT edited). Same filters/shrinkage/significance rule, but the
# design and the study subset are arguments, so a single study can be fitted with ~ group.
# usage: deseq2_design.R <counts.tsv> <sample_map.tsv (srr,label,study,group)> <gene_id_to_symbol.tsv> <prefix> <min_count> <min_samples> <padj> <lfc> <design: group|study+group> <study|ALL> [exclude_srr,comma]
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(jsonlite) })
a <- commandArgs(trailingOnly = TRUE)
prefix <- a[4]; min_count <- as.numeric(a[5]); min_samples <- as.numeric(a[6]); padj_cut <- as.numeric(a[7]); lfc_cut <- as.numeric(a[8]); design_s <- a[9]; stud <- a[10]
excl <- if (length(a) >= 11 && nzchar(a[11])) strsplit(a[11], ",")[[1]] else character(0)
ss <- read.delim(a[2], stringsAsFactors = FALSE); rownames(ss) <- ss$srr
symbols <- read.delim(a[3], row.names = 1, stringsAsFactors = FALSE)
cn <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
keep_s <- intersect(colnames(cn), rownames(ss)); keep_s <- setdiff(keep_s, excl)
if (stud != "ALL") keep_s <- keep_s[ss[keep_s, "study"] == stud]
cn <- cn[, keep_s, drop = FALSE]; ss <- ss[keep_s, , drop = FALSE]
ss$study <- factor(ss$study); ss$group <- relevel(factor(ss$group), ref = "normal")
print(table(ss$study, ss$group))
if (nlevels(ss$group) != 2) stop("ERROR: need both normal and tumour samples")
des <- switch(design_s, "group" = ~ group, "study+group" = ~ study + group, stop("unknown design"))
dds <- DESeqDataSetFromMatrix(round(cn), colData = ss, design = des)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
message(sprintf("pre-filter: >= %g counts in >= %g samples: kept %d of %d genes", min_count, min_samples, sum(keep), length(keep)))
dds <- dds[keep, ]; dds <- DESeq(dds)
coef <- "group_tumour_vs_normal"; stopifnot(coef %in% resultsNames(dds))
res <- results(dds, name = coef); shr <- lfcShrink(dds, coef = coef, type = "apeglm", res = res)
out <- as.data.frame(res); out$log2FoldChange_shrunk <- shr$log2FoldChange; out$gene_id <- rownames(out)
out$gene_name <- symbols[out$gene_id, "gene_name"]; out$gene_type <- symbols[out$gene_id, "gene_type"]
out <- out[order(out$padj), c("gene_id", "gene_name", "gene_type", "baseMean", "log2FoldChange", "log2FoldChange_shrunk", "lfcSE", "stat", "pvalue", "padj")]
w <- function(d, f) write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
w(out, paste0(prefix, "_results_all.tsv"))
sig <- !is.na(out$padj) & out$padj < padj_cut & abs(out$log2FoldChange) >= lfc_cut
w(out[sig & out$log2FoldChange > 0, ], paste0(prefix, "_significant_up.tsv")); w(out[sig & out$log2FoldChange < 0, ], paste0(prefix, "_significant_down.tsv"))
w(data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE), paste0(prefix, "_normalized_counts.tsv"))
write_json(list(design = design_s, study = stud, reference = "normal", samples = rownames(ss), n_samples = ncol(dds), n_genes_after_prefilter = nrow(dds),
  n_tested = sum(!is.na(out$padj)), n_padj_lt_cut = sum(!is.na(out$padj) & out$padj < padj_cut), padj_cut = padj_cut, lfc_cut = lfc_cut,
  n_sig_up = sum(sig & out$log2FoldChange > 0), n_sig_down = sum(sig & out$log2FoldChange < 0)), paste0(prefix, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), paste0(prefix, "_sessionInfo.txt"))
