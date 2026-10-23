#!/usr/bin/env Rscript
# Verify part 2: Step 1 (S2 read counts), Step 3 (pathways vs S4/S5/S6), Step 4b (leave-one-out T2_B = SRR13401202).
suppressPackageStartupMessages({ library(jsonlite); library(DESeq2); library(apeglm) })
B <- "/workspace/Biomarker_pipeline"; R <- file.path(B, "results"); V <- file.path(R, "09_verification"); O <- file.path(V, "step9_paper_comparison"); X <- file.path(V, "input/MOESM2_sheets")
ss <- read.delim(file.path(R, "sample_map.tsv"), stringsAsFactors = FALSE); rownames(ss) <- ss$srr
# ---------- Step 1
s2 <- readLines(file.path(X, "Table_S2.tsv")); s2 <- strsplit(s2, "\t"); s2 <- s2[vapply(s2, function(x) grepl("^SRR[0-9]+$", x[1]), NA)]
S2 <- data.frame(srr = vapply(s2, `[`, "", 1), s2_reads = as.numeric(gsub("[^0-9].*$", "", gsub(",", "", vapply(s2, `[`, "", 2)))), stringsAsFactors = FALSE)
t1 <- ss[, c("srr", "label", "study")]
t1$our_fastp_pairs <- vapply(t1$srr, function(s) fromJSON(file.path(R, "01_fastp", s, "fastp.json"))$summary$before_filtering$total_reads / 2, 0)
t1$s2_reads <- S2$s2_reads[match(t1$srr, S2$srr)]; t1$diff <- t1$our_fastp_pairs - t1$s2_reads; t1$match <- ifelse(is.na(t1$s2_reads), "not in S2", ifelse(t1$diff == 0, "yes", "no"))
write.table(t1, file.path(O, "step1_read_counts_vs_S2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(sprintf("runs compared: %d; match yes: %d; differ: %d; not in S2: %d", nrow(t1), sum(t1$match == "yes"), sum(t1$match == "no"), sum(t1$match == "not in S2")),
  capture.output(print(t1[, c("srr", "label", "our_fastp_pairs", "s2_reads", "match")], row.names = FALSE, format = "d", big.mark = ","))), file.path(O, "step1_summary.txt"))
# ---------- Step 3
rd <- function(f) read.delim(file.path(X, f), skip = 1, stringsAsFactors = FALSE, quote = "")
stripgo <- function(x) trimws(sub("\\s*\\(GO:[0-9]+\\)\\s*$", "", x))
norm <- function(x) tolower(gsub("[^a-z0-9]+", " ", tolower(trimws(x))))
cmp <- function(ours, paper, tag) {
  on <- unique(norm(ours$name)); pn <- unique(norm(paper$name)); both <- intersect(on, pn)
  d <- list(both = both, paper_only = setdiff(pn, on), ours_only = setdiff(on, pn))
  show <- function(n, df, ord) { k <- df$name[match(n, norm(df$name))]; k[order(match(n, norm(df$name[order(df[[ord]])])))][1:min(10, length(n))] }
  ours$k <- norm(ours$name); paper$k <- norm(paper$name)
  top <- function(keys, df) { x <- df[df$k %in% keys, ]; x <- x[order(x$p), ]; head(unique(x$name), 10) }
  c(sprintf("== %s: ours significant %d terms, paper significant %d terms; in both %d; paper only %d; ours only %d", tag, length(on), length(pn), length(both), length(d$paper_only), length(d$ours_only)),
    "  top 10 in both (by our padj): ", paste0("    ", top(both, ours)), "  top 10 paper only (by paper padj): ", paste0("    ", top(d$paper_only, paper)), "  top 10 ours only (by our padj): ", paste0("    ", top(d$ours_only, ours)), "")
}
out <- c()
for (dir in c("up", "down")) {
  pa <- rd(if (dir == "up") "Table_S6_Upregulated_.tsv" else "Table_S6_Downregulated_.tsv"); pa <- pa[pa$Adjusted.P.value < 0.05, ]
  ou <- read.delim(file.path(R, "07_pathway/ora", sprintf("ora_%s_KEGG.tsv", dir)), stringsAsFactors = FALSE); ou <- ou[ou$p.adjust < 0.05, ]
  out <- c(out, cmp(data.frame(name = ou$Description, p = ou$p.adjust), data.frame(name = pa$Term, p = pa$Adjusted.P.value), sprintf("KEGG %s", dir)))
  pg <- rd(if (dir == "up") "Table_S5_Biological_process_.tsv" else "Table_S4_Biological_process_.tsv"); pg <- pg[pg$Adjusted.P.value < 0.05, ]
  og <- read.delim(file.path(R, "07_pathway/ora", sprintf("ora_%s_GO_BP.tsv", dir)), stringsAsFactors = FALSE); og <- og[og$p.adjust < 0.05, ]
  out <- c(out, cmp(data.frame(name = og$Description, p = og$p.adjust), data.frame(name = stripgo(pg$Term), p = pg$Adjusted.P.value), sprintf("GO BP %s", dir)))
}
writeLines(c("padj < 0.05 in both; terms matched by name (lower case, punctuation removed). Paper: Enrichr-style tables on its integrated 4 SRA studies + TCGA READ gene lists; ours: clusterProfiler ORA (KEGG REST release at run time, GO BP) on our 2-study DESeq2 up/down sets. KEGG/GO releases and background gene sets differ.", "", out), file.path(O, "step3_pathway_comparison.txt"))
# ---------- Step 4b leave-one-out
pub <- read.delim(file.path(R, "06_dge/deseq2_results_all.tsv"), stringsAsFactors = FALSE); rownames(pub) <- pub$gene_id
fc <- as.matrix(read.delim(file.path(R, "04_counts/counts_featurecounts.tsv"), row.names = 1, check.names = FALSE))
keep <- setdiff(colnames(fc), "SRR13401202"); cd <- ss[keep, ]; cd$study <- factor(cd$study); cd$group <- relevel(factor(cd$group), "normal")
d <- DESeqDataSetFromMatrix(round(fc[, keep]), cd, ~ study + group); d <- d[rowSums(counts(d) >= 10) >= 5, ]; d <- DESeq(d, quiet = TRUE)
r <- results(d, name = "group_tumour_vs_normal"); a <- as.data.frame(r); a$gene_id <- rownames(a); a$gene_name <- pub[a$gene_id, "gene_name"]
a$log2FoldChange_shrunk <- lfcShrink(d, coef = "group_tumour_vs_normal", type = "apeglm", res = r, quiet = TRUE)$log2FoldChange
write.table(a, file.path(O, "step4b_deseq2_without_T2_B_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sig <- function(o) !is.na(o$padj) & o$padj < 0.05 & abs(o$log2FoldChange) >= 1
sa <- a$gene_id[sig(a)]; sm0 <- pub$gene_id[sig(pub)]; ov <- intersect(sa, sm0); cm <- intersect(a$gene_id, pub$gene_id); rownames(a) <- a$gene_id
S3 <- read.delim(file.path(V, "input/hosseini_2023_TableS3_genes.tsv"), stringsAsFactors = FALSE); S3$gene <- trimws(S3$gene); S3 <- S3[!duplicated(S3$gene), ]
ag <- function(o, only_sig) { o <- o[!is.na(o$log2FoldChange) & o$log2FoldChange != 0 & (!only_sig | sig(o)), ]; x <- merge(o, S3, by.x = "gene_name", by.y = "gene"); sprintf("%d genes, %.1f%%", nrow(x), 100 * mean(sign(x$log2FoldChange) == ifelse(x$direction_in_cancer == "up", 1, -1))) }
g18 <- c("GUCA2A","CLCA4","CLCA1","AQP8","SLC26A3","APOA4","PYY","GUCA2B","ZG16","CD177","COL3A1","COL1A1","SPARC","TGFBI","CXCL8","KRT6A","KRT6B","KRT17")
t18 <- data.frame(gene = g18, main_lfc = pub$log2FoldChange[match(g18, pub$gene_name)], main_padj = pub$padj[match(g18, pub$gene_name)], wo_lfc = a$log2FoldChange[match(g18, a$gene_name)], wo_padj = a$padj[match(g18, a$gene_name)]); t18[, -1] <- signif(t18[, -1], 3)
q <- c("Leave-one-out: SRR13401202 (T2_B) removed, 18 samples, design ~ study + group, ref normal, apeglm, prefilter >= 10 counts in >= 5 samples, featureCounts route", capture.output(print(table(cd$study, cd$group))),
  sprintf("tested genes %d (main %d)", sum(!is.na(a$padj)), sum(!is.na(pub$padj))),
  sprintf("significant (padj<0.05, |lfc|>=1 unshrunk): without T2_B %d; main %d; overlap %d; Jaccard %.3f; same direction %d of %d", length(sa), length(sm0), length(ov), length(ov) / length(union(sa, sm0)), sum(sign(a[ov, "log2FoldChange"]) == sign(pub[ov, "log2FoldChange"])), length(ov)),
  sprintf("log2FC over %d common genes: Spearman %.3f, Pearson %.3f", length(cm), cor(a[cm, "log2FoldChange"], pub[cm, "log2FoldChange"], method = "spearman", use = "complete.obs"), cor(a[cm, "log2FoldChange"], pub[cm, "log2FoldChange"], use = "complete.obs")),
  sprintf("Table S3 direction, all S3-listed genes tested: without T2_B %s; main %s", ag(a, FALSE), ag(pub, FALSE)),
  sprintf("Table S3 direction, our significant genes: without T2_B %s; main %s", ag(a, TRUE), ag(pub, TRUE)), "", capture.output(print(t18, row.names = FALSE)))
writeLines(q, file.path(O, "step4b_leave_one_out_T2_B.txt"))
