#!/usr/bin/env Rscript
# Steps 2, 3, 5, 6 (DESeq2 based). Refits the pipeline model (same prefilter, design ~ study + group, reference normal) to get the dds object
# the pipeline did not save, checks it reproduces results/06_dge, then: diagnostics (2), per-study replication (3), Hosseini S3 (5), sex check (6).
# Cutoffs are the pipeline's (padj 0.05, |log2FC| >= 1, prefilter >= 10 counts in >= 5 samples); nothing is tuned.
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(ggplot2) })
B <- "/workspace/Biomarker_pipeline"; R <- file.path(B, "results"); V <- file.path(R, "09_verification")
PADJ <- 0.05; LFC <- 1; MINC <- 10; MINS <- 5
ss <- read.delim(file.path(R, "sample_map.tsv"), stringsAsFactors = FALSE); rownames(ss) <- ss$srr
sym <- read.delim(file.path(B, "ref/gene_id_to_symbol.tsv"), row.names = 1, stringsAsFactors = FALSE)
fc <- as.matrix(read.delim(file.path(R, "04_counts/counts_featurecounts.tsv"), row.names = 1, check.names = FALSE))
txi <- readRDS(file.path(R, "04_counts/salmon_tximport.rds"))
pub_fc <- read.delim(file.path(R, "06_dge/deseq2_results_all.tsv"), stringsAsFactors = FALSE)
pub_sal <- read.delim(file.path(R, "06_dge/deseq2_salmon_results_all.tsv"), stringsAsFactors = FALSE)
`%+%` <- function(a, b) paste0(a, b)
sink(file.path(V, "step2_deseq2_diagnostics/log_fits.txt"), split = TRUE)

mk <- function(cn, txi = NULL, samples, design, min_samples = MINS) {
  cd <- ss[samples, , drop = FALSE]; cd$study <- factor(cd$study); cd$group <- relevel(factor(cd$group), ref = "normal")
  dds <- if (!is.null(txi)) { t2 <- txi; for (n in c("counts", "abundance", "length")) t2[[n]] <- txi[[n]][, samples, drop = FALSE]; DESeqDataSetFromTximport(t2, colData = cd, design = design) }
         else DESeqDataSetFromMatrix(round(cn[, samples, drop = FALSE]), colData = cd, design = design)
  dds <- dds[rowSums(counts(dds) >= MINC) >= min_samples, ]; DESeq(dds, quiet = TRUE)
}
tab <- function(dds, shrink = TRUE) { res <- results(dds, name = "group_tumour_vs_normal"); o <- as.data.frame(res); o$gene_id <- rownames(o)
  o$gene_name <- sym[o$gene_id, "gene_name"]
  if (shrink) o$lfc_shrunk <- lfcShrink(dds, coef = "group_tumour_vs_normal", type = "apeglm", res = res, quiet = TRUE)$log2FoldChange
  list(tab = o, res = res) }
sigf <- function(o) !is.na(o$padj) & o$padj < PADJ & abs(o$log2FoldChange) >= LFC

# ---------------------------------------------------------------- Step 2
D2 <- file.path(V, "step2_deseq2_diagnostics"); summ <- c()
for (route in c("featureCounts", "Salmon")) {
  samples <- if (route == "Salmon") colnames(txi$counts) else colnames(fc)
  dds <- if (route == "Salmon") mk(NULL, txi, samples, ~ study + group) else mk(fc, NULL, samples, ~ study + group)
  x <- tab(dds); o <- x$tab; res <- x$res; pub <- if (route == "Salmon") pub_sal else pub_fc; rownames(pub) <- pub$gene_id
  common <- intersect(o$gene_id, pub$gene_id)
  repro <- c(n_refit = nrow(o), n_pipeline = nrow(pub), max_abs_diff_log2FC = max(abs(o[common, "log2FoldChange"] - pub[common, "log2FoldChange"]), na.rm = TRUE),
             max_abs_diff_shrunk = max(abs(o[common, "lfc_shrunk"] - pub[common, "log2FoldChange_shrunk"]), na.rm = TRUE),
             padj_same_NA_pattern = all(is.na(o[common, "padj"]) == is.na(pub[common, "padj"])),
             max_rel_diff_padj = max(abs(o[common, "padj"] - pub[common, "padj"]) / pub[common, "padj"], na.rm = TRUE))
  tag <- gsub("featureCounts", "fc", route)
  png(file.path(D2, sprintf("%s_pvalue_histogram.png", tag)), 900, 650, res = 120)
  hist(o$pvalue[!is.na(o$pvalue)], breaks = 50, col = "grey70", border = "white", main = sprintf("%s: p-value histogram (tumour vs normal, ~ study + group)", route), xlab = "p-value"); dev.off()
  png(file.path(D2, sprintf("%s_dispersion.png", tag)), 900, 700, res = 120); plotDispEsts(dds, main = sprintf("%s: dispersion estimates", route)); dev.off()
  d <- o[!is.na(o$log2FoldChange) & !is.na(o$lfc_shrunk), ]
  ggsave(file.path(D2, sprintf("%s_shrunk_vs_unshrunk.png", tag)), ggplot(d, aes(log2FoldChange, lfc_shrunk, colour = !is.na(padj) & padj < PADJ)) + geom_point(size = 0.5, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) + scale_colour_manual(values = c("grey60", "red3"), name = "padj < 0.05") + theme_bw() +
    labs(title = sprintf("%s: apeglm-shrunk vs unshrunk log2FC", route), x = "unshrunk log2FC", y = "apeglm log2FC"), width = 6.5, height = 6, dpi = 130)
  na_p <- sum(is.na(o$pvalue)); na_padj <- sum(is.na(o$padj)); indep <- sum(is.na(o$padj) & !is.na(o$pvalue))
  cooks_out <- sum(is.na(o$pvalue) & o$baseMean > 0)
  s <- sigf(o); sh <- sum(s & abs(o$lfc_shrunk) < LFC, na.rm = TRUE)
  summ <- c(summ, sprintf("== %s (%d samples)", route, ncol(dds)),
    sprintf("reproduces pipeline results/06_dge: genes %d vs %d; max |diff| log2FC %.2e, shrunk %.2e; padj NA pattern identical: %s; max relative padj difference %.2e",
            repro["n_refit"], repro["n_pipeline"], repro["max_abs_diff_log2FC"], repro["max_abs_diff_shrunk"], repro["padj_same_NA_pattern"], repro["max_rel_diff_padj"]),
    sprintf("genes after prefilter: %d; padj NA: %d (%.1f%%): pvalue NA from Cook's-distance outlier flag (baseMean > 0): %d; independent filtering (pvalue present, padj NA): %d; independent-filter baseMean threshold: %.2f",
            nrow(o), na_padj, 100 * na_padj / nrow(o), cooks_out, indep, metadata(res)$filterThreshold),
    sprintf("p-value histogram: fraction of p < 0.05 = %.3f, fraction of p > 0.95 = %.3f (flat = 0.05 in each bin of width 0.05)", mean(o$pvalue < 0.05, na.rm = TRUE), mean(o$pvalue > 0.95, na.rm = TRUE)),
    sprintf("dispersion: median gene-wise %.3f, median fitted %.3f, median final %.3f; fraction of genes with dispersion outlier flag: %.3f",
            median(mcols(dds)$dispGeneEst, na.rm = TRUE), median(mcols(dds)$dispFit, na.rm = TRUE), median(dispersions(dds), na.rm = TRUE), mean(mcols(dds)$dispOutlier, na.rm = TRUE)),
    sprintf("shrinkage: Pearson r(unshrunk, apeglm) = %.4f over %d genes; significant genes (padj < 0.05 & |unshrunk| >= 1) whose shrunk |log2FC| < 1: %d of %d; median |shrunk|/|unshrunk| among significant: %.3f",
            cor(d$log2FoldChange, d$lfc_shrunk), nrow(d), sh, sum(s, na.rm = TRUE), median(abs(o$lfc_shrunk[s]) / abs(o$log2FoldChange[s]), na.rm = TRUE)), "")
  if (route == "featureCounts") { dds_fc <- dds; res_fc <- o } else { res_sal <- o }
}
writeLines(summ, file.path(D2, "deseq2_diagnostics_summary.txt")); cat(summ, sep = "\n")

# ---------------------------------------------------------------- Step 3: replication per study (featureCounts route)
D3 <- file.path(V, "step3_replication")
genes18 <- c("GUCA2A", "CLCA4", "CLCA1", "AQP8", "SLC26A3", "APOA4", "PYY", "GUCA2B", "ZG16", "CD177", "COL3A1", "COL1A1", "SPARC", "TGFBI", "CXCL8", "KRT6A", "KRT6B", "KRT17")
exp18 <- c(rep("down", 10), rep("up", 8))
stud <- function(st) { sm <- rownames(ss)[ss$study == st & rownames(ss) %in% colnames(fc)]; cat(sprintf("study %s: %d tumour, %d normal\n", st, sum(ss[sm, "group"] == "tumour"), sum(ss[sm, "group"] == "normal")))
  tab(mk(fc, NULL, sm, ~ group), shrink = FALSE)$tab }
A <- stud("SRP219837"); Bs <- stud("SRP301216"); F <- res_fc
rownames(A) <- A$gene_id; rownames(Bs) <- Bs$gene_id; rownames(F) <- F$gene_id
cm <- intersect(A$gene_id[!is.na(A$log2FoldChange)], Bs$gene_id[!is.na(Bs$log2FoldChange)])
r_all <- cor(A[cm, "log2FoldChange"], Bs[cm, "log2FoldChange"]); rs_all <- cor(A[cm, "log2FoldChange"], Bs[cm, "log2FoldChange"], method = "spearman")
sa <- sigf(A); sb <- sigf(Bs); both <- intersect(A$gene_id[sa], Bs$gene_id[sb]); same <- both[sign(A[both, "log2FoldChange"]) == sign(Bs[both, "log2FoldChange"])]
either <- union(A$gene_id[sa], Bs$gene_id[sb]); r_either <- cor(A[either, "log2FoldChange"], Bs[either, "log2FoldChange"], use = "complete.obs")
padj_both <- cm[!is.na(A[cm, "padj"]) & !is.na(Bs[cm, "padj"])]
sigA_dir <- A[padj_both, "padj"] < PADJ; sigB_dir <- Bs[padj_both, "padj"] < PADJ
sameA_B <- sum(sigA_dir & sigB_dir & sign(A[padj_both, "log2FoldChange"]) == sign(Bs[padj_both, "log2FoldChange"]))
g18 <- do.call(rbind, lapply(seq_along(genes18), function(i) { g <- genes18[i]; pick <- function(o) { z <- o[!is.na(o$gene_name) & o$gene_name == g, ]; if (!nrow(z)) c(NA, NA) else { z <- z[which.max(z$baseMean), ]; c(z$log2FoldChange, z$padj) } }
  a <- pick(A); b <- pick(Bs); f <- pick(F); data.frame(gene = g, expected = exp18[i], log2FC_studyA = a[1], padj_studyA = a[2], log2FC_studyB = b[1], padj_studyB = b[2], log2FC_both_studies = f[1], padj_both_studies = f[2]) }))
dirm <- function(l, e) ifelse(is.na(l), NA, (l > 0) == (e == "up"))
g18$dir_match_A <- dirm(g18$log2FC_studyA, g18$expected); g18$dir_match_B <- dirm(g18$log2FC_studyB, g18$expected); g18$dir_match_both <- dirm(g18$log2FC_both_studies, g18$expected)
g18$sig_A <- !is.na(g18$padj_studyA) & g18$padj_studyA < PADJ; g18$sig_B <- !is.na(g18$padj_studyB) & g18$padj_studyB < PADJ
write.table(g18, file.path(D3, "hosseini18_per_study.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_id = cm, gene_name = A[cm, "gene_name"], log2FC_A = A[cm, "log2FoldChange"], padj_A = A[cm, "padj"], log2FC_B = Bs[cm, "log2FoldChange"], padj_B = Bs[cm, "padj"]), file.path(D3, "per_study_log2FC_common_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ggsave(file.path(D3, "log2FC_studyA_vs_studyB.png"), ggplot(data.frame(a = A[cm, "log2FoldChange"], b = Bs[cm, "log2FoldChange"], s = ifelse(A[cm, "gene_id"] %in% both, "sig in both", "other")), aes(a, b, colour = s)) + geom_point(size = 0.5, alpha = 0.5) +
  scale_colour_manual(values = c("sig in both" = "red3", other = "grey60")) + geom_abline(slope = 1, intercept = 0, linetype = 2) + theme_bw() + labs(x = "log2FC study A (5 T vs 5 N)", y = "log2FC study B (4 T vs 5 N)", title = sprintf("per-study log2FC, r = %.3f (%d genes)", r_all, length(cm))), width = 6.5, height = 6, dpi = 130)
s3 <- c(sprintf("study A (SRP219837, 5 tumour vs 5 normal): genes tested (padj not NA) %d, significant (padj < %g & |log2FC| >= %g) %d (up %d, down %d)", sum(!is.na(A$padj)), PADJ, LFC, sum(sa), sum(sa & A$log2FoldChange > 0), sum(sa & A$log2FoldChange < 0)),
        sprintf("study B (SRP301216, 4 tumour vs 5 normal): genes tested %d, significant %d (up %d, down %d)", sum(!is.na(Bs$padj)), sum(sb), sum(sb & Bs$log2FoldChange > 0), sum(sb & Bs$log2FoldChange < 0)),
        sprintf("log2FC correlation A vs B over %d common genes: Pearson %.4f, Spearman %.4f; over genes significant in either study (%d): Pearson %.4f", length(cm), r_all, rs_all, length(either), r_either),
        sprintf("genes significant in BOTH studies (each with padj < %g & |log2FC| >= %g): %d; same direction: %d (%.1f%%); union %d, Jaccard %.3f", PADJ, LFC, length(both), length(same), 100 * length(same) / max(1, length(both)), length(either), length(both) / length(either)),
        sprintf("looser view, padj < %g in both (no log2FC cutoff), same direction: %d of %d genes tested in both", PADJ, sameA_B, sum(sigA_dir & sigB_dir)),
        sprintf("overlap of the per-study significant genes with the full 19-sample significant set (%d genes): study A %d, study B %d", sum(sigf(F)), sum(A$gene_id[sa] %in% F$gene_id[sigf(F)]), sum(Bs$gene_id[sb] %in% F$gene_id[sigf(F)])),
        sprintf("18 Hosseini genes: direction matches expected in study A %d/18, study B %d/18, both studies together %d/18; significant (padj < %g): study A %d, study B %d", sum(g18$dir_match_A, na.rm = TRUE), sum(g18$dir_match_B, na.rm = TRUE), sum(g18$dir_match_both, na.rm = TRUE), PADJ, sum(g18$sig_A), sum(g18$sig_B)),
        "", capture.output(print(g18[, c("gene", "expected", "log2FC_studyA", "padj_studyA", "log2FC_studyB", "padj_studyB", "dir_match_A", "dir_match_B")], row.names = FALSE, digits = 3)))
writeLines(s3, file.path(D3, "replication_summary.txt")); cat(s3, sep = "\n")

# ---------------------------------------------------------------- Step 5: Hosseini Table S3 (directional only)
D5 <- file.path(V, "step5_hosseini_S3")
S3 <- read.delim(file.path(V, "input/hosseini_2023_TableS3_genes.tsv"), stringsAsFactors = FALSE)
S3$gene <- trimws(S3$gene); S3 <- S3[!duplicated(S3$gene), ]
O <- F[!is.na(F$gene_name) & F$gene_name != "" & !is.na(F$padj), ]; O <- O[order(-O$baseMean), ]; O <- O[!duplicated(O$gene_name), ]    # one row per symbol (highest baseMean)
N <- nrow(O); O$ours <- ifelse(O$log2FoldChange > 0, "up", "down"); O$ours_sig <- ifelse(sigf(O), O$ours, "ns")
m <- merge(O[, c("gene_name", "baseMean", "log2FoldChange", "padj", "ours", "ours_sig")], S3, by.x = "gene_name", by.y = "gene")
m$agree <- m$ours == m$direction_in_cancer
write.table(m, file.path(D5, "S3_genes_tested_in_our_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
hyp <- function(a, b, label) { k <- length(intersect(a, b)); K <- length(b); n <- length(a)
  c(test = label, ours = n, paper_in_universe = K, overlap = k, expected = round(n * K / N, 1), fold = round(k / (n * K / N), 2), p_enrichment = signif(phyper(k - 1, K, N - K, n, lower.tail = FALSE), 3)) }
su <- O$gene_name[O$ours_sig == "up"]; sd <- O$gene_name[O$ours_sig == "down"]; pu <- m$gene_name[m$direction_in_cancer == "up"]; pd <- m$gene_name[m$direction_in_cancer == "down"]
H <- as.data.frame(rbind(hyp(su, pu, "ours UP & padj/lfc  vs  S3 up"), hyp(sd, pd, "ours DOWN & padj/lfc  vs  S3 down"), hyp(su, pd, "ours UP vs S3 down (discordant)"), hyp(sd, pu, "ours DOWN vs S3 up (discordant)")), stringsAsFactors = FALSE)
write.table(H, file.path(D5, "hypergeometric_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dir_tab <- function(d, lab) { x <- m[m$direction_in_cancer == d, ]; c(sprintf("S3 %s genes in our tested set: %d; same direction (sign of our log2FC): %d (%.1f%%)", lab, nrow(x), sum(x$agree), 100 * mean(x$agree)),
  { y <- x[x$ours_sig != "ns", ]; sprintf("   of those, significant in ours (padj < %g, |log2FC| >= %g): %d; same direction %d (%.1f%%)", PADJ, LFC, nrow(y), sum(y$agree), 100 * mean(y$agree)) }) }
bt <- function(d) { x <- m[m$direction_in_cancer == d, ]; signif(binom.test(sum(x$agree), nrow(x), p = if (d == "up") mean(O$log2FoldChange > 0) else mean(O$log2FoldChange < 0))$p.value, 3) }
s5 <- c(sprintf("Hosseini & Nemati 2023, Supplementary Table S3: %d up, %d down genes listed (full DEG list, tumour vs normal).", sum(S3$direction_in_cancer == "up"), sum(S3$direction_in_cancer == "down")),
        sprintf("Our tested genes with a symbol: %d (genes with padj not NA, one row per symbol). S3 genes found among them: up %d of %d, down %d of %d.", N, length(pu), sum(S3$direction_in_cancer == "up"), length(pd), sum(S3$direction_in_cancer == "down")),
        dir_tab("up", "UP"), dir_tab("down", "DOWN"),
        sprintf("binomial test of direction agreement against the share of our tested genes that go that way (up %.1f%% / down %.1f%% of our genes): S3 up p = %s, S3 down p = %s", 100 * mean(O$log2FoldChange > 0), 100 * mean(O$log2FoldChange < 0), bt("up"), bt("down")),
        sprintf("overall direction agreement for all S3 genes tested: %d of %d (%.1f%%)", sum(m$agree), nrow(m), 100 * mean(m$agree)),
        "", "Hypergeometric overlap (universe = our tested genes with symbol, N = " %+% N %+% "; paper set restricted to the universe):", capture.output(print(H, row.names = FALSE)),
        "", "Caveats: S3 lists only the paper's DEGs (4,752 up / 285 down: very asymmetric), so the test asks whether our significant genes are over-represented among them; the paper's own tested universe is unknown. Directional only: no cutoff changed.")
writeLines(s5, file.path(D5, "S3_agreement_summary.txt")); cat(s5, sep = "\n")

# ---------------------------------------------------------------- Step 6: sex check
D6 <- file.path(V, "step6_sex_check")
nc <- read.delim(file.path(R, "06_dge/deseq2_normalized_counts.tsv"), row.names = 1, check.names = FALSE)
ids <- sapply(c("XIST", "RPS4Y1", "DDX3Y", "KDM5D"), function(g) { z <- F[!is.na(F$gene_name) & F$gene_name == g, ]; z$gene_id[which.max(z$baseMean)] })
L <- log2(t(nc[ids, , drop = FALSE]) + 1); colnames(L) <- names(ids)
L <- as.data.frame(L); L$srr <- rownames(L); L$label <- ss[L$srr, "label"]; L$study <- ss[L$srr, "study"]; L$group <- ss[L$srr, "group"]
L$Y_mean <- rowMeans(L[, c("RPS4Y1", "DDX3Y", "KDM5D")])
L$call <- ifelse(L$Y_mean >= 3, "male-like", ifelse(L$XIST >= 5, "female-like", "ambiguous"))
# patient identity: study A from the ENA sample titles (COAD-T004 / COAD-N004 = patient 004; label numbers T1_A, N1_A are NOT patient numbers);
# study B from the sample titles too ([PT1] / [PN1]); the label number equals the patient number there.
tA <- read.delim(file.path(V, "input/studyA_SRP219837_ena_sample_titles.tsv"), stringsAsFactors = FALSE); rownames(tA) <- tA$srr
L$patient <- ifelse(L$study == "SRP219837", sub(".*COAD-[TN]([0-9]+).*", "A\\1", tA[L$srr, "sample_title"]), paste0("B", sub("^[TN]([0-9]).*", "\\1", L$label)))
pat <- aggregate(call ~ patient, L, function(x) paste(sort(unique(x)), collapse = "/")); pat$n_samples <- as.vector(table(L$patient)[pat$patient])
pat$consistent <- !grepl("/", pat$call)
L <- L[order(L$study, L$patient, L$group), c("srr", "label", "study", "group", "patient", "XIST", "RPS4Y1", "DDX3Y", "KDM5D", "Y_mean", "call")]
write.table(L, file.path(D6, "sex_markers_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pat, file.path(D6, "sex_call_per_patient_pair.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
Lp <- L; Lp$label <- factor(Lp$label, levels = Lp$label)
ggsave(file.path(D6, "sex_markers.png"), ggplot(Lp, aes(XIST, Y_mean, colour = study, shape = group, label = label)) + geom_point(size = 3) + geom_text(vjust = -1, size = 3, show.legend = FALSE) + theme_bw() +
  geom_hline(yintercept = 3, linetype = 2, colour = "grey50") + labs(x = "XIST  log2(normalized count + 1)", y = "mean log2(norm + 1) of RPS4Y1, DDX3Y, KDM5D", title = "sex-marker expression (19 samples)"), width = 7.5, height = 6, dpi = 130)
s6 <- c("log2(DESeq2 normalized count + 1) of XIST (female) and RPS4Y1/DDX3Y/KDM5D (male, non-PAR Y genes; PAR genes were removed from the gene set):", capture.output(print(L, row.names = FALSE, digits = 3)), "",
        "call rule (descriptive, no metadata used): male-like if mean Y-gene value >= 3; female-like if < 3 and XIST >= 5; else ambiguous", "calls: " %+% "", capture.output(print(table(L$study, L$call))), "",
        "pairs (tumour and adjacent normal of the same patient number within a study should have the same sex):", capture.output(print(pat, row.names = FALSE)))
writeLines(s6, file.path(D6, "sex_check_summary.txt")); cat(s6, sep = "\n")
sink()
