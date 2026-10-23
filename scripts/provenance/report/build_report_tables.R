#!/usr/bin/env Rscript
# Builds the summary tables and the ROC plot of diagnostic_report_20261015/ from EXISTING files only (read-only on inputs).
# No model is refitted; no DESeq2, LASSO or pathway analysis is run. 19 development samples only.
# Run from the project root: env/prefix/r_classifier/bin/Rscript results_rc1prep_19/diagnostic_report_20261015/build_report_tables.R
suppressPackageStartupMessages(library(pROC))
R <- "results_rc1prep_19"; O <- file.path(R, "diagnostic_report_20261015")
wt <- function(x, f) write.table(x, file.path(O, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
log_lines <- character(0); say <- function(...) { m <- paste0(...); message(m); log_lines <<- c(log_lines, m) }

# ---- 1. DE summary per route (from the DESeq2 result tables; rule 1 = padj < 0.05 and |shrunk log2FC| >= 1)
de_route <- function(f, route) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  r1 <- !is.na(d$padj) & d$padj < 0.05 & !is.na(d$log2FoldChange_shrunk) & abs(d$log2FoldChange_shrunk) >= 1
  un <- !is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1
  data.frame(route = route, source_file = f, genes_in_table_after_prefilter = nrow(d),
             genes_tested_nonNA_pvalue = sum(!is.na(d$pvalue)), genes_nonNA_padj = sum(!is.na(d$padj)),
             padj_lt_0.05 = sum(!is.na(d$padj) & d$padj < 0.05),
             rule1_shrunk_total = sum(r1), rule1_shrunk_up = sum(r1 & d$log2FoldChange_shrunk > 0), rule1_shrunk_down = sum(r1 & d$log2FoldChange_shrunk < 0),
             unshrunk_lfc1_total = sum(un), unshrunk_lfc1_up = sum(un & d$log2FoldChange > 0), unshrunk_lfc1_down = sum(un & d$log2FoldChange < 0))
}
de <- rbind(de_route(file.path(R, "deseq2_results_all.tsv"), "featureCounts (primary)"),
            de_route(file.path(R, "deseq2_salmon_results_all.tsv"), "Salmon/tximport"))
wt(de, "de_summary_by_route.tsv")

# ---- 2. Candidate funnel (18 of 19) and final 1,390 split
ct <- read.delim(file.path(R, "rule3b_rerun18_cfg/candidate_table.tsv"), stringsAsFactors = FALSE)
fin <- read.delim(file.path(R, "final_candidates/final_candidates_18of19.tsv"), stringsAsFactors = FALSE)
istrue <- function(x) x %in% c(TRUE, "TRUE")
pass_ids <- ct$gene_id[istrue(ct$all_rules_pass)]
say("final_candidates_18of19.tsv rows: ", nrow(fin), "; candidate table (rule3b_rerun18_cfg) all_rules_pass: ", length(pass_ids),
    "; identical gene_id sets: ", setequal(pass_ids, fin$gene_id))
pc <- fin$gene_biotype == "protein_coding"
fs <- data.frame(item = c("final candidates", "up (direction_fc)", "down (direction_fc)", "protein_coding", "protein_coding up", "protein_coding down",
                          "baseMean_ge_50_both TRUE", "loo_n_pass = 19", "loo_n_pass = 18"),
                 n = c(nrow(fin), sum(fin$direction_fc == "up"), sum(fin$direction_fc == "down"), sum(pc), sum(pc & fin$direction_fc == "up"),
                       sum(pc & fin$direction_fc == "down"), sum(istrue(fin$baseMean_ge_50_both)), sum(fin$loo_n_pass == 19), sum(fin$loo_n_pass == 18)))
say("direction_fc values in final list: ", paste(names(table(fin$direction_fc)), table(fin$direction_fc), collapse = "; "))
wt(fs, "final_candidates_split.tsv")
file.copy(file.path(R, "rule3b_rerun18_cfg/funnel_counts.tsv"), file.path(O, "candidate_funnel_18of19.tsv"))
file.copy(file.path(R, "final_candidates/final_candidates_18of19.tsv"), file.path(O, "final_candidates_1390_full.tsv"))

# ---- 3. Top candidates (featureCounts route): by padj, and by |shrunk log2FC|
cols <- c("gene_id", "symbol", "gene_biotype", "direction_fc", "baseMean_fc", "log2FC_shrunk_fc", "padj_fc", "log2FC_shrunk_salmon", "padj_salmon", "loo_n_pass")
by_p <- fin[order(fin$padj_fc, -abs(fin$log2FC_shrunk_fc)), cols][1:25, ]
by_l <- fin[order(-abs(fin$log2FC_shrunk_fc), fin$padj_fc), cols][1:25, ]
wt(cbind(rank = 1:25, by_p), "top25_candidates_by_padj.tsv")
wt(cbind(rank = 1:25, by_l), "top25_candidates_by_abs_shrunk_log2FC.tsv")

# ---- 4. Reactome ORA / GSEA top 10 (values copied from the result tables, not recomputed)
ora_top <- function(f, lab) {
  d <- read.delim(file.path(R, "reactome", f), stringsAsFactors = FALSE); d <- d[order(d$p.adjust, d$pvalue), ]
  h <- head(d, 10); data.frame(set = lab, rank = seq_len(nrow(h)), pathway = h$ID, GeneRatio = h$GeneRatio, BgRatio = h$BgRatio,
                               p.adjust_BH = h$p.adjust, qvalue = h$qvalue, Count = h$Count, significant_BH_lt_0.05 = h$significant_BH_lt_0.05)
}
ora <- rbind(ora_top("ora_up.tsv", "up"), ora_top("ora_down.tsv", "down"), ora_top("ora_combined.tsv", "combined"))
wt(ora, "reactome_ora_top10.tsv")
g <- read.delim(file.path(R, "reactome/gsea_reactome.tsv"), stringsAsFactors = FALSE)
gt <- function(x, lab) { x <- x[order(x$padj, x$pval), ]; h <- head(x, 10)
  data.frame(direction = lab, rank = seq_len(nrow(h)), pathway = h$pathway, NES = h$NES, padj = h$padj, size = h$size, significant_padj_lt_0.05 = h$significant_padj_lt_0.05) }
gs <- rbind(gt(g[g$NES > 0, ], "positive NES (higher in tumour)"), gt(g[g$NES < 0, ], "negative NES (lower in tumour)"))
wt(gs, "reactome_gsea_top10.tsv")
sig <- istrue(g$significant_padj_lt_0.05)
say("GSEA significant: ", sum(sig), " (NES>0: ", sum(sig & g$NES > 0), ", NES<0: ", sum(sig & g$NES < 0), ")")
for (f in c("analysis_summary.tsv", "mapping_summary.tsv", "tool_versions.txt")) file.copy(file.path(R, "reactome", f), file.path(O, paste0("reactome_", f)))

# ---- 5. Locked model: copy files; recompute development scores from the locked genes/formula (no refit) and compare with the existing score file
lk <- "lock/diagnostic_v1"
file.copy(file.path(lk, "final_model_genes.tsv"), file.path(O, "locked_model_genes.tsv"))
mg <- read.delim(file.path(lk, "final_model_genes.tsv"), stringsAsFactors = FALSE)
CUT <- 1.352325902   # final_model_spec.txt line "cutoff = 1.352325902"
spec <- readLines(file.path(lk, "final_model_spec.txt")); stopifnot(any(grepl("cutoff = 1.352325902", spec, fixed = TRUE)))
cn <- as.matrix(read.delim("results/04_counts/counts_featurecounts.tsv", row.names = 1, check.names = FALSE)); stopifnot(ncol(cn) == 19)
stopifnot(all(mg$gene_id %in% rownames(cn)))
L <- log2(sweep(cn, 2, colSums(cn), "/") * 1e6 + 1)
up <- mg$gene_id[as.numeric(mg$sign) > 0]; dn <- mg$gene_id[as.numeric(mg$sign) < 0]
sc <- colMeans(L[up, , drop = FALSE]) - colMeans(L[dn, , drop = FALSE])
ex <- read.delim(file.path(R, "lock_proposal_v2/development_scores.tsv"), stringsAsFactors = FALSE)
ex$score_recomputed_locked_model <- sc[ex$sample]; ex$abs_diff <- abs(ex$score - ex$score_recomputed_locked_model)
ex$pred_tumour_recomputed <- as.integer(ex$score_recomputed_locked_model > CUT)
say("locked-model score recompute: n_up ", length(up), ", n_down ", length(dn), "; max |existing - recomputed| = ", signif(max(ex$abs_diff), 3),
    "; predictions identical: ", all(ex$pred_tumour == ex$pred_tumour_recomputed))
wt(ex, "development_scores_locked_model.tsv")
y <- as.integer(ex$truth == "tumour")
roc1 <- roc(y, ex$score, levels = c(0, 1), direction = "<", quiet = TRUE); ci1 <- ci.auc(roc1, method = "delong")
tp <- sum(ex$pred_tumour == 1 & y == 1); tn <- sum(ex$pred_tumour == 0 & y == 0)
perf <- data.frame(n = nrow(ex), n_tumour = sum(y), n_normal = sum(1 - y), AUC = as.numeric(auc(roc1)), AUC_ci_low_delong = ci1[1], AUC_ci_high_delong = ci1[3],
                   TP = tp, FN = sum(y) - tp, TN = tn, FP = sum(1 - y) - tn, cutoff = CUT,
                   min_tumour_score = min(ex$score[y == 1]), max_normal_score = max(ex$score[y == 0]),
                   label = "development data, optimistic, not validation")
wt(perf, "development_performance_locked_model.tsv")
png(file.path(O, "roc_locked_model_development.png"), width = 1600, height = 1600, res = 250)
par(mar = c(5, 5, 6, 2))
plot(1 - roc1$specificities, roc1$sensitivities, type = "s", lwd = 2.5, col = "#1f5fa8", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - specificity", ylab = "sensitivity", main = "")
abline(0, 1, lty = 2, col = "grey60")
mtext("Locked 9-gene score, lock-diagnostic-v1 (n = 19)", side = 3, line = 3.6, font = 2, cex = 0.95)
mtext("DEVELOPMENT DATA, OPTIMISTIC, NOT VALIDATION", side = 3, line = 2.2, col = "#b00020", font = 2, cex = 0.95)
mtext(sprintf("AUC %.3f (CI not informative at AUC 1); genes + cutoff chosen on these samples", perf$AUC), side = 3, line = 0.8, cex = 0.75)
dev.off()
writeLines(log_lines, file.path(O, "build_report_tables.log"))
