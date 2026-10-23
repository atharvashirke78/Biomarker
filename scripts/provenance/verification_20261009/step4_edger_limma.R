#!/usr/bin/env Rscript
# Step 4: second method on the same featureCounts matrix, ~ study + group, same prefilter genes as DESeq2 (>= 10 counts in >= 5 samples).
# edgeR (TMM, glm quasi-likelihood F-test) and limma-voom. Significance as in the pipeline: FDR < 0.05 and |log2FC| >= 1. Nothing tuned.
suppressPackageStartupMessages({ library(edgeR); library(limma) })
B <- "/workspace/Biomarker_pipeline"; R <- file.path(B, "results"); V <- file.path(R, "09_verification"); D <- file.path(V, "step4_edger_limma")
ss <- read.delim(file.path(R, "sample_map.tsv"), stringsAsFactors = FALSE); rownames(ss) <- ss$srr
fc <- as.matrix(read.delim(file.path(R, "04_counts/counts_featurecounts.tsv"), row.names = 1, check.names = FALSE))
ss <- ss[colnames(fc), ]; ss$study <- factor(ss$study); ss$group <- relevel(factor(ss$group), ref = "normal")
sym <- read.delim(file.path(B, "ref/gene_id_to_symbol.tsv"), row.names = 1, stringsAsFactors = FALSE)
ds <- read.delim(file.path(R, "06_dge/deseq2_results_all.tsv"), stringsAsFactors = FALSE); rownames(ds) <- ds$gene_id
keep <- rowSums(round(fc) >= 10) >= 5; cn <- round(fc[keep, ]); stopifnot(setequal(rownames(cn), ds$gene_id))
design <- model.matrix(~ study + group, ss); co <- "grouptumour"
y <- DGEList(cn); y <- calcNormFactors(y); y <- estimateDisp(y, design)
fit <- glmQLFit(y, design); qlf <- glmQLFTest(fit, coef = co); et <- topTags(qlf, n = Inf, sort.by = "none")$table
v <- voom(y, design, plot = FALSE); lf <- eBayes(lmFit(v, design), robust = TRUE); lt <- topTable(lf, coef = co, number = Inf, sort.by = "none")
ed <- data.frame(gene_id = rownames(et), gene_name = sym[rownames(et), "gene_name"], edgeR_log2FC = et$logFC, edgeR_FDR = et$FDR, edgeR_PValue = et$PValue, stringsAsFactors = FALSE)
ed$limma_log2FC <- lt[ed$gene_id, "logFC"]; ed$limma_FDR <- lt[ed$gene_id, "adj.P.Val"]; ed$limma_t <- lt[ed$gene_id, "t"]
ed$deseq2_log2FC <- ds[ed$gene_id, "log2FoldChange"]; ed$deseq2_padj <- ds[ed$gene_id, "padj"]; ed$deseq2_stat <- ds[ed$gene_id, "stat"]
write.table(ed, file.path(D, "edgeR_limma_vs_deseq2_all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sg <- function(l, p) !is.na(p) & p < 0.05 & abs(l) >= 1
S <- list(DESeq2 = sg(ed$deseq2_log2FC, ed$deseq2_padj), edgeR = sg(ed$edgeR_log2FC, ed$edgeR_FDR), limma_voom = sg(ed$limma_log2FC, ed$limma_FDR))
jac <- function(a, b) sum(a & b) / sum(a | b)
cmp <- function(x, lx, y, ly) { ok <- !is.na(lx) & !is.na(ly); both <- x & y; dir <- sum(both & sign(lx) == sign(ly))
  c(n_x = sum(x), n_y = sum(y), overlap = sum(both), Jaccard = round(jac(x, y), 3), same_dir_of_overlap = dir, pearson_log2FC_all = round(cor(lx[ok], ly[ok]), 4), spearman_log2FC_all = round(cor(lx[ok], ly[ok], method = "spearman"), 4),
    pearson_log2FC_sig_either = round(cor(lx[x | y], ly[x | y], use = "complete.obs"), 4)) }
T <- rbind("DESeq2 vs edgeR" = cmp(S$DESeq2, ed$deseq2_log2FC, S$edgeR, ed$edgeR_log2FC), "DESeq2 vs limma-voom" = cmp(S$DESeq2, ed$deseq2_log2FC, S$limma_voom, ed$limma_log2FC),
           "edgeR vs limma-voom" = cmp(S$edgeR, ed$edgeR_log2FC, S$limma_voom, ed$limma_log2FC))
write.table(cbind(comparison = rownames(T), T), file.path(D, "method_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ct <- c(sprintf("19 samples, %d genes (same prefilter as DESeq2), design ~ study + group, reference normal; significance FDR/padj < 0.05 & |log2FC| >= 1", nrow(ed)),
        sprintf("significant: DESeq2 %d, edgeR %d, limma-voom %d; padj/FDR < 0.05 only (no log2FC cutoff): DESeq2 %d, edgeR %d, limma-voom %d", sum(S$DESeq2), sum(S$edgeR), sum(S$limma_voom),
                sum(ed$deseq2_padj < 0.05, na.rm = TRUE), sum(ed$edgeR_FDR < 0.05), sum(ed$limma_FDR < 0.05)),
        sprintf("genes significant in all three: %d", sum(S$DESeq2 & S$edgeR & S$limma_voom)), "", capture.output(print(T)), "",
        sprintf("edgeR common dispersion %.3f (BCV %.3f); TMM norm factors range %.2f-%.2f", y$common.dispersion, sqrt(y$common.dispersion), min(y$samples$norm.factors), max(y$samples$norm.factors)),
        sprintf("DESeq2 padj NA (Cook's outlier flag) %d genes: edgeR FDR < 0.05 for %d of them", sum(is.na(ed$deseq2_padj)), sum(is.na(ed$deseq2_padj) & ed$edgeR_FDR < 0.05)))
writeLines(ct, file.path(D, "method_comparison_summary.txt")); cat(ct, sep = "\n")
sc <- function(f, yv, lab, r) { png(file.path(D, f), 800, 800, res = 130); plot(ed$deseq2_log2FC, yv, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3), xlab = "DESeq2 log2FC (unshrunk)", ylab = lab, main = sprintf("r = %.4f", r)); abline(0, 1, col = "red3", lty = 2); dev.off() }
sc("log2FC_deseq2_vs_edgeR.png", ed$edgeR_log2FC, "edgeR logFC", T[1, "pearson_log2FC_all"]); sc("log2FC_deseq2_vs_limma_voom.png", ed$limma_log2FC, "limma-voom logFC", T[2, "pearson_log2FC_all"])
