#!/usr/bin/env Rscript
# Extra questions: (2) outlier genes, (3) refit without N4_A and N5_A. Same prefilter/design/cutoffs as the pipeline; reads only.
suppressPackageStartupMessages({ library(DESeq2); library(apeglm) })
B <- "/workspace/Biomarker_pipeline"; R <- file.path(B, "results"); O <- file.path(R, "09_verification/step8_extra_questions")
ss <- read.delim(file.path(R, "sample_map.tsv"), stringsAsFactors = FALSE); rownames(ss) <- ss$srr
fc <- as.matrix(read.delim(file.path(R, "04_counts/counts_featurecounts.tsv"), row.names = 1, check.names = FALSE))
pub <- read.delim(file.path(R, "06_dge/deseq2_results_all.tsv"), stringsAsFactors = FALSE)
nc <- read.delim(file.path(R, "06_dge/deseq2_normalized_counts.tsv"), row.names = 1, check.names = FALSE)
genes <- c("KRT6A", "KRT6B", "KRT17", "CXCL8")
# ---- Q2
q2 <- pub[pub$gene_name %in% genes, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "log2FoldChange_shrunk", "padj")]
q2$padj_NA <- is.na(q2$padj)
sm <- ss[colnames(nc), ]
lines <- capture.output(print(q2, row.names = FALSE))
for (g in q2$gene_id) { v <- unlist(nc[g, ]); names(v) <- paste0(sm[names(v), "label"], "(", substr(sm[names(v), "group"], 1, 1), ")")
  lines <- c(lines, "", paste("normalized counts (top 5):", pub$gene_name[pub$gene_id == g]), capture.output(print(round(sort(v, decreasing = TRUE)[1:5], 1))))}
# Cook's flag check: refit full model, look at assays cooks / replaceable
cn <- fc[, rownames(ss)[rownames(ss) %in% colnames(fc)]]
fit <- function(samples) { cd <- ss[samples, ]; cd$study <- factor(cd$study); cd$group <- relevel(factor(cd$group), "normal")
  d <- DESeqDataSetFromMatrix(round(fc[, samples]), cd, ~ study + group); d <- d[rowSums(counts(d) >= 10) >= 5, ]; d <- DESeq(d, quiet = TRUE)
  r <- results(d, name = "group_tumour_vs_normal"); o <- as.data.frame(r); o$gene_id <- rownames(o); o$gene_name <- pub$gene_name[match(o$gene_id, pub$gene_id)]
  list(d = d, o = o) }
full <- fit(colnames(fc))
cook <- assays(full$d)[["cooks"]]
for (g in q2$gene_id) { cv <- sort(cook[g, ], decreasing = TRUE)[1:2]
  lines <- c(lines, sprintf("%s: max Cook's %.2f in %s; pvalue %s, padj %s -> %s", pub$gene_name[pub$gene_id == g], cv[1], sm[names(cv)[1], "label"],
   format(full$o[g, "pvalue"], digits = 3), format(full$o[g, "padj"], digits = 3), if (is.na(full$o[g, "pvalue"])) "pvalue NA (all-zero/outlier)" else if (is.na(full$o[g, "padj"])) "pvalue present, padj NA (Cook's outlier flag or independent filtering)" else "padj present")) }
writeLines(lines, file.path(O, "q2_outlier_genes.txt"))
# ---- Q3
drop <- ss$srr[ss$label %in% c("N4_A", "N5_A")]; keep <- setdiff(colnames(fc), drop)
sub <- fit(keep); a <- sub$o; rownames(a) <- a$gene_id; m <- pub; rownames(m) <- m$gene_id
sig <- function(o) !is.na(o$padj) & o$padj < 0.05 & abs(o$log2FoldChange) >= 1
sa <- rownames(a)[sig(a)]; sm0 <- rownames(m)[sig(m)]; ov <- intersect(sa, sm0); cm <- intersect(rownames(a), rownames(m))
S3 <- read.delim(file.path(R, "09_verification/input/hosseini_2023_TableS3_genes.tsv"), stringsAsFactors = FALSE); S3$gene <- trimws(S3$gene); S3 <- S3[!duplicated(S3$gene), ]
agree <- function(o) { x <- merge(o[!is.na(o$log2FoldChange) & o$log2FoldChange != 0, ], S3, by.x = "gene_name", by.y = "gene")
  c(n = nrow(x), pct = 100 * mean(sign(x$log2FoldChange) == ifelse(x$direction_in_cancer == "up", 1, -1)))}
ag <- function(o) { x <- merge(o[sig(o), ], S3, by.x = "gene_name", by.y = "gene"); c(n = nrow(x), pct = 100 * mean(sign(x$log2FoldChange) == ifelse(x$direction_in_cancer == "up", 1, -1))) }
q3 <- c(sprintf("dropped: %s", paste(ss[drop, "label"], collapse = ", ")), sprintf("samples %d (main 19); tested genes %d (main %d)", length(keep), sum(!is.na(a$padj)), sum(!is.na(m$padj))),
 table(ss[keep, "study"], ss[keep, "group"]) |> capture.output(),
 sprintf("significant (padj<0.05 & |lfc|>=1, unshrunk): without N4_A+N5_A %d; main %d; overlap %d; Jaccard %.3f", length(sa), length(sm0), length(ov), length(ov) / length(union(sa, sm0))),
 sprintf("same direction among overlap: %d of %d", sum(sign(a[ov, "log2FoldChange"]) == sign(m[ov, "log2FoldChange"])), length(ov)),
 sprintf("log2FC over %d common genes: Spearman %.3f, Pearson %.3f", length(cm), cor(a[cm, "log2FoldChange"], m[cm, "log2FoldChange"], method = "spearman", use = "complete.obs"), cor(a[cm, "log2FoldChange"], m[cm, "log2FoldChange"], use = "complete.obs")),
 sprintf("Table S3 direction agreement, all S3-listed genes tested: without %d genes %.1f%%; main %d genes %.1f%%", agree(a)["n"], agree(a)["pct"], agree(m)["n"], agree(m)["pct"]),
 sprintf("Table S3 direction agreement among our significant genes: without %d genes %.1f%%; main %d genes %.1f%%", ag(a)["n"], ag(a)["pct"], ag(m)["n"], ag(m)["pct"]))
g18 <- c("GUCA2A","CLCA4","CLCA1","AQP8","SLC26A3","APOA4","PYY","GUCA2B","ZG16","CD177","COL3A1","COL1A1","SPARC","TGFBI","CXCL8","KRT6A","KRT6B","KRT17")
t18 <- data.frame(gene = g18, main_lfc = m$log2FoldChange[match(g18, m$gene_name)], main_padj = m$padj[match(g18, m$gene_name)], wo_lfc = a$log2FoldChange[match(g18, a$gene_name)], wo_padj = a$padj[match(g18, a$gene_name)])
t18[, -1] <- signif(t18[, -1], 3)
writeLines(c(q3, "", capture.output(print(t18, row.names = FALSE))), file.path(O, "q3_without_N4A_N5A.txt"))
write.table(a, file.path(O, "deseq2_without_N4A_N5A_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
