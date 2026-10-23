#!/usr/bin/env Rscript
# Step 7: paired sensitivity analysis. Study A: the 4 real tumour/normal pairs (patients 004, 012, 014, 017 from the ENA sample titles), design ~ patient + group,
# reference normal, featureCounts route, same prefilter and significance rules as the pipeline. Extra (labelled): study B paired (patients 1-4 + lone normal 5).
# Compared with (a) the main unpaired 19-sample result in results/06_dge and (b) the unpaired study-A-only fit (5 vs 5) from step 3.
suppressPackageStartupMessages(library(DESeq2))
B <- "/workspace/Biomarker_pipeline"; R <- file.path(B, "results"); V <- file.path(R, "09_verification"); D <- file.path(V, "step7_paired_sensitivity")
pm <- read.delim(file.path(V, "sample_map_with_patient.tsv"), stringsAsFactors = FALSE); rownames(pm) <- pm$srr
fc <- as.matrix(read.delim(file.path(R, "04_counts/counts_featurecounts.tsv"), row.names = 1, check.names = FALSE))
main <- read.delim(file.path(R, "06_dge/deseq2_results_all.tsv"), stringsAsFactors = FALSE); rownames(main) <- main$gene_id
sym <- read.delim(file.path(B, "ref/gene_id_to_symbol.tsv"), row.names = 1, stringsAsFactors = FALSE)
sig <- function(l, p) !is.na(p) & p < 0.05 & abs(l) >= 1
fit <- function(samples) { cd <- pm[samples, c("srr", "study", "group", "patient_id")]; cd$group <- relevel(factor(cd$group), "normal"); cd$patient <- factor(cd$patient_id)
  dds <- DESeqDataSetFromMatrix(round(fc[, samples]), cd, ~ patient + group); dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]; dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "group_tumour_vs_normal"); o <- as.data.frame(res); o$gene_id <- rownames(o); o$gene_name <- sym[o$gene_id, "gene_name"]; o }
two <- names(which(tapply(pm$group, pm$patient_id, function(x) length(unique(x))) == 2))
sA <- pm$srr[pm$study == "SRP219837" & pm$patient_id %in% two]; sB <- pm$srr[pm$study == "SRP301216" & pm$patient_id %in% two & pm$srr %in% colnames(fc)]
sB <- c(sB, pm$srr[pm$study == "SRP301216" & pm$srr %in% colnames(fc) & !(pm$srr %in% sB)])  # lone normal N5_B stays (absorbed by its own patient level, no effect on the group estimate)
cat("study A paired samples:", paste(pm[sA, "label"], collapse = " "), "\n"); cat("study B samples:", paste(pm[sB, "label"], collapse = " "), "\n")
PA <- fit(sA); PB <- fit(sB)
UA <- {  # unpaired study A only, 5 vs 5, ~ group (as in step 3)
  sa <- pm$srr[pm$study == "SRP219837"]; cd <- pm[sa, ]; cd$group <- relevel(factor(cd$group), "normal"); d <- DESeqDataSetFromMatrix(round(fc[, sa]), cd, ~ group); d <- d[rowSums(counts(d) >= 10) >= 5, ]; d <- DESeq(d, quiet = TRUE)
  o <- as.data.frame(results(d, name = "group_tumour_vs_normal")); o$gene_id <- rownames(o); o$gene_name <- sym[o$gene_id, "gene_name"]; o }
write.table(PA, file.path(D, "studyA_paired_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); write.table(PB, file.path(D, "studyB_paired_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cmp <- function(P, Q, nm) { cm <- intersect(P$gene_id[!is.na(P$log2FoldChange)], Q$gene_id[!is.na(Q$log2FoldChange)]); sp <- sig(P$log2FoldChange, P$padj); names(sp) <- P$gene_id; sq <- sig(Q$log2FoldChange, Q$padj); names(sq) <- Q$gene_id
  a <- names(sp)[sp]; b <- names(sq)[sq]; both <- intersect(a, b); same <- both[sign(P[match(both, P$gene_id), "log2FoldChange"]) == sign(Q[match(both, Q$gene_id), "log2FoldChange"])]
  Pl <- P$log2FoldChange[match(cm, P$gene_id)]; Ql <- Q$log2FoldChange[match(cm, Q$gene_id)]
  data.frame(comparison = nm, genes_compared = length(cm), n_sig_paired = length(a), n_sig_other = length(b), overlap = length(both), same_direction = length(same), Jaccard = round(length(both) / length(union(a, b)), 3),
             Spearman_log2FC = round(cor(Pl, Ql, method = "spearman"), 3), Pearson_log2FC = round(cor(Pl, Ql), 3)) }
T <- rbind(cmp(PA, main, "study A paired (4 pairs) vs main unpaired 19-sample"), cmp(PA, UA, "study A paired (4 pairs) vs study A unpaired (5 v 5)"), cmp(PB, main, "[extra] study B paired vs main unpaired 19-sample"))
write.table(T, file.path(D, "paired_vs_unpaired_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
genes <- c("GUCA2A","CLCA4","CLCA1","AQP8","SLC26A3","APOA4","PYY","GUCA2B","ZG16","CD177","COL3A1","COL1A1","SPARC","TGFBI","CXCL8","KRT6A","KRT6B","KRT17"); ex <- c(rep("down", 10), rep("up", 8))
pick <- function(o, g) { z <- o[!is.na(o$gene_name) & o$gene_name == g, ]; if (!nrow(z)) c(NA, NA) else { z <- z[which.max(z$baseMean), ]; c(z$log2FoldChange, z$padj) } }
G <- do.call(rbind, lapply(seq_along(genes), function(i) { a <- pick(PA, genes[i]); u <- pick(UA, genes[i]); m <- pick(main, genes[i]); b <- pick(PB, genes[i])
  data.frame(gene = genes[i], expected = ex[i], log2FC_A_paired = a[1], padj_A_paired = a[2], log2FC_A_unpaired = u[1], padj_A_unpaired = u[2], log2FC_main_19 = m[1], padj_main_19 = m[2], log2FC_B_paired = b[1], padj_B_paired = b[2]) }))
G$dir_A_paired <- ifelse(is.na(G$log2FC_A_paired), NA, (G$log2FC_A_paired > 0) == (G$expected == "up")); G$sig_A_paired <- !is.na(G$padj_A_paired) & G$padj_A_paired < 0.05
write.table(G, file.path(D, "hosseini18_paired_vs_unpaired.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 3))
txt <- c(sprintf("Study A paired: %d samples (%s); design ~ patient + group, ref normal; prefilter >= 10 counts in >= 5 samples: %d genes tested (padj not NA %d)", length(sA), paste(pm[sA, "label"], collapse = ", "), nrow(PA), sum(!is.na(PA$padj))),
         sprintf("Study B (extra): %d samples; %d genes", length(sB), nrow(PB)),
         sprintf("significant (padj < 0.05 & |log2FC| >= 1): A paired %d (up %d, down %d); A unpaired %d; main 19-sample %d; B paired %d", sum(sig(PA$log2FoldChange, PA$padj)), sum(sig(PA$log2FoldChange, PA$padj) & PA$log2FoldChange > 0), sum(sig(PA$log2FoldChange, PA$padj) & PA$log2FoldChange < 0),
                 sum(sig(UA$log2FoldChange, UA$padj)), sum(sig(main$log2FoldChange, main$padj)), sum(sig(PB$log2FoldChange, PB$padj))), "",
         capture.output(print(T, row.names = FALSE)), "", "18 Hosseini genes (log2FC / padj):",
         capture.output(print(data.frame(gene = G$gene, exp = G$expected, A_paired = paste0(fmt(G$log2FC_A_paired), " / ", fmt(G$padj_A_paired)), A_unpaired = paste0(fmt(G$log2FC_A_unpaired), " / ", fmt(G$padj_A_unpaired)),
                                         main_19 = paste0(fmt(G$log2FC_main_19), " / ", fmt(G$padj_main_19)), B_paired = paste0(fmt(G$log2FC_B_paired), " / ", fmt(G$padj_B_paired))), row.names = FALSE, right = FALSE)), "",
         sprintf("18 genes, direction as expected in study A paired: %d of %d testable (significant padj < 0.05: %d)", sum(G$dir_A_paired, na.rm = TRUE), sum(!is.na(G$dir_A_paired)), sum(G$sig_A_paired)))
writeLines(txt, file.path(D, "paired_sensitivity_summary.txt")); cat(txt, sep = "\n")
