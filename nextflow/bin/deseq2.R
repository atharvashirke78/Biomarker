#!/usr/bin/env Rscript
# REPLACED (copy of /workspace/biomarker/pipeline/bin/deseq2.R was ~ condition, tximport only).
# usage: deseq2.R <matrix|tximport> <counts.tsv|.rds> <sample_sheet.tsv (sample_id,label,group,baseline,study,batch,count_source)> <gene_id_to_symbol.tsv> <prefix> <min_count> <min_samples (number, or "auto" = smallest group in the sample sheet)> <alpha (padj cut; also results() alpha)> <lfc> <n_volcano_labels> <n_heatmap_genes> [<volcano_label_seed>]; sample IDs = SRR
# 2026-10-19: optional 12th argument = seed of the volcano label placement (conf/input.config de_volcano_seed; owner decision); absent = unseeded as before.
# Refactor (Step 2): group, baseline group and batch are read from the sample sheet (no hard-coded names); design ~ batch + group (~ group when there is a single batch).
# Input checks (baseline, group sizes, integers, sheet/matrix match, estimable design) are done before this script by validate_analysis_input.py.
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(ggplot2); library(ggrepel); library(pheatmap); library(jsonlite) })
a <- commandArgs(trailingOnly = TRUE)
mode <- a[1]; prefix <- a[5]; min_count <- as.numeric(a[6]); alpha <- as.numeric(a[8]); padj_cut <- alpha; lfc_cut <- as.numeric(a[9]); n_labels <- as.integer(a[10]); n_heat <- as.integer(a[11]); vseed <- if (length(a) >= 12) as.integer(a[12]) else NA
ss <- read.delim(a[3], stringsAsFactors = FALSE, colClasses = "character"); rownames(ss) <- ss$sample_id
symbols <- read.delim(a[4], row.names = 1, stringsAsFactors = FALSE)
if (mode == "tximport") { txi <- readRDS(a[2]); cn <- txi$counts } else { cn <- as.matrix(read.delim(a[2], row.names = 1, check.names = FALSE)) }
ss <- ss[colnames(cn), , drop = FALSE]
base_group <- unique(ss$group[toupper(ss$baseline) == "TRUE"]); stopifnot(length(base_group) == 1)
ss$study <- factor(ss$study); ss$batch <- factor(ss$batch); ss$group <- relevel(factor(ss$group), ref = base_group)
print(table(ss$batch, ss$group))
if (nlevels(ss$group) != 2) stop("ERROR: need exactly two groups")
if (any(table(ss$batch, ss$group) == 0)) stop("ERROR: a batch x group cell is empty; ~ batch + group is not estimable")
design <- if (nlevels(ss$batch) > 1) ~ batch + group else ~ group
design_txt <- if (nlevels(ss$batch) > 1) "~ batch + group" else "~ group"
dds <- if (mode == "tximport") DESeqDataSetFromTximport(txi, colData = ss, design = design) else DESeqDataSetFromMatrix(round(cn), colData = ss, design = design)
min_samples <- if (a[7] == "auto") min(table(ss$group)) else as.numeric(a[7])  # smallest group size from the sample sheet (DESeq2 vignette pre-filter)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
message(sprintf("pre-filter: >= %g counts in >= %g samples: kept %d of %d genes", min_count, min_samples, sum(keep), length(keep)))
dds <- dds[keep, ]; dds <- DESeq(dds)
test_group <- setdiff(levels(ss$group), base_group)
coef <- paste0("group_", test_group, "_vs_", base_group); stopifnot(coef %in% resultsNames(dds))
res <- results(dds, name = coef, alpha = alpha); shr <- lfcShrink(dds, coef = coef, type = "apeglm", res = res)
out <- as.data.frame(res); out$log2FoldChange_shrunk <- shr$log2FoldChange; out$gene_id <- rownames(out)
out$gene_name <- symbols[out$gene_id, "gene_name"]; out$gene_type <- symbols[out$gene_id, "gene_type"]
out <- out[order(out$padj), c("gene_id", "gene_name", "gene_type", "baseMean", "log2FoldChange", "log2FoldChange_shrunk", "lfcSE", "stat", "pvalue", "padj")]
w <- function(d, f) write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
w(out, paste0(prefix, "_results_all.tsv"))
sig <- !is.na(out$padj) & out$padj < padj_cut & abs(out$log2FoldChange) >= lfc_cut
w(out[sig & out$log2FoldChange > 0, ], paste0(prefix, "_significant_up.tsv")); w(out[sig & out$log2FoldChange < 0, ], paste0(prefix, "_significant_down.tsv"))
w(data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE), paste0(prefix, "_normalized_counts.tsv"))
v <- out[!is.na(out$padj), ]
v$cls <- ifelse(v$padj < padj_cut & v$log2FoldChange >= lfc_cut, "up", ifelse(v$padj < padj_cut & v$log2FoldChange <= -lfc_cut, "down", "ns"))
lab <- head(v[v$cls != "ns", ], n_labels)
ggsave(paste0(prefix, "_volcano.png"), ggplot(v, aes(log2FoldChange, -log10(padj), colour = cls)) + geom_point(size = 0.6, alpha = 0.6) +
  scale_colour_manual(values = c(up = "red3", down = "blue3", ns = "grey60")) +
  geom_text_repel(data = lab, aes(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)), size = 3, colour = "black", max.overlaps = 30, seed = vseed) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2) + geom_hline(yintercept = -log10(padj_cut), linetype = 2) + theme_bw() +
  labs(title = paste0(prefix, ": ", test_group, " vs ", base_group, " (", design_txt, ")")), width = 7, height = 6, dpi = 150)
ggsave(paste0(prefix, "_MA.png"), ggplot(data.frame(b = out$baseMean, l = out$log2FoldChange_shrunk, s = !is.na(out$padj) & out$padj < padj_cut), aes(b, l, colour = s)) +
  geom_point(size = 0.6, alpha = 0.6) + scale_x_log10() + scale_colour_manual(values = c("grey60", "red3")) + geom_hline(yintercept = 0) + theme_bw() + labs(y = "apeglm log2FC"), width = 6, height = 5, dpi = 150)
vsd <- tryCatch(vst(dds, blind = FALSE), error = function(e) varianceStabilizingTransformation(dds, blind = FALSE))
top <- head(out$gene_id[!is.na(out$padj)], n_heat); mm <- assay(vsd)[top, , drop = FALSE]; mm <- mm - rowMeans(mm)
rn <- symbols[top, "gene_name"]; rownames(mm) <- make.unique(ifelse(is.na(rn) | rn == "", top, rn))
pheatmap(mm, annotation_col = ss[, c("batch", "group")], filename = paste0(prefix, "_top50_heatmap.png"), width = 8, height = 10, fontsize_row = 7)
write_json(list(route = mode, design = design_txt, reference = base_group, tested_group = test_group, n_samples = ncol(dds), n_genes_after_prefilter = nrow(dds),
  n_tested = sum(!is.na(out$padj)), n_padj_lt_cut = sum(!is.na(out$padj) & out$padj < padj_cut), padj_cut = padj_cut, lfc_cut = lfc_cut,
  n_sig_up = sum(sig & out$log2FoldChange > 0), n_sig_down = sum(sig & out$log2FoldChange < 0)), paste0(prefix, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), paste0(prefix, "_sessionInfo.txt"))
