#!/usr/bin/env Rscript
# DEVELOPMENT-ONLY relative-expression analysis on the 19 discovery samples. NOT part of the locked diagnostic model
# (lock-diagnostic-v1 normalises by featureCounts library total only; this analysis must not be substituted into it).
#
# Expression measure:  E = log2(DESeq2 median-of-ratios normalised count + PC), PC = 1 (primary)
#   normalised counts are REUSED from the 19-sample DESeq2 fit (results_rc1prep_19/deseq2_normalized_counts.tsv)
# Relative expression: RE(sample, gene) = E(gene) - mean over the 10 frozen reference genes (config/reference_genes.yaml) of E(ref)
#   An RNA-seq log-ratio, NOT a qPCR delta-Ct (no cycle threshold, no amplification efficiency). Gene length differs between
#   genes, so RE is comparable between samples for the same gene only (the tumour-vs-normal contrast cancels length).
# Statistics per gene: lm(RE ~ study + group) -> tumour - normal (log2), 95% CI, p; BH within the 9 model genes and, separately,
#   within the 1,390 candidates. Sensitivity: paired within-patient analysis (8 tumour/normal pairs; mean difference, t CI,
#   exact Wilcoxon signed-rank), pseudocounts 0.5 and 4, and raw counts normalised by the reference genes only.
# usage (project root): env/prefix/r_dge/bin/Rscript scripts/downstream/relative_expression_dev.R <new out dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: relative_expression_dev.R <new out dir>")
out <- args[1]
if (file.exists(out)) stop("ERROR: output folder exists (never overwritten): ", out)
R <- "results_rc1prep_19"

# ---- inputs (read-only)
nc <- as.matrix(read.delim(file.path(R, "deseq2_normalized_counts.tsv"), row.names = 1, check.names = FALSE))
raw <- as.matrix(read.delim("results/04_counts/counts_featurecounts.tsv", row.names = 1, check.names = FALSE))[rownames(nc), colnames(nc)]
ss <- read.delim(file.path(R, "sample_sheet_rc1prep.tsv"), stringsAsFactors = FALSE); ss <- ss[match(colnames(nc), ss$sample_id), ]
stopifnot(identical(ss$sample_id, colnames(nc)), ncol(nc) == 19)
pm <- read.delim("results/09_verification/sample_map_with_patient.tsv", stringsAsFactors = FALSE)
ss$patient <- pm$patient_id[match(ss$sample_id, pm$srr)]; stopifnot(!anyNA(ss$patient))
ss$group <- relevel(factor(ss$group), ref = "normal"); ss$study <- factor(ss$study)
ref_ids <- sub("^.*gene_id:\\s*(ENSG[0-9]+).*$", "\\1", grep("gene_id: ENSG", readLines("config/reference_genes.yaml"), value = TRUE))
unv <- sub("\\..*$", "", rownames(nc))
ref_rows <- rownames(nc)[match(ref_ids, unv)]; stopifnot(length(ref_ids) == 10, !anyNA(ref_rows))
model <- read.delim("lock/diagnostic_v1/final_model_genes.tsv", stringsAsFactors = FALSE); stopifnot(all(model$gene_id %in% rownames(nc)))
cand <- read.delim(file.path(R, "final_candidates/final_candidates_18of19.tsv"), stringsAsFactors = FALSE); stopifnot(all(cand$gene_id %in% rownames(nc)))
de <- read.delim(file.path(R, "deseq2_results_all.tsv"), stringsAsFactors = FALSE)
sym <- read.delim("ref/gene_id_to_symbol.tsv", stringsAsFactors = FALSE)
ref_sym <- sym$gene_name[match(ref_rows, sym$gene_id)]

# ---- relative expression
re_of <- function(mat, pc) { E <- log2(mat + pc); sweep(E, 2, colMeans(E[ref_rows, , drop = FALSE]), "-") }
RE <- re_of(nc, 1)
refmean <- colMeans(log2(nc[ref_rows, ] + 1))

fit_gene <- function(y) {
  f <- lm(y ~ ss$study + ss$group); co <- summary(f)$coefficients; ci <- confint(f)
  k <- "ss$grouptumour"; c(est = co[k, 1], se = co[k, 2], ci_low = ci[k, 1], ci_high = ci[k, 2], p = co[k, 4],
                          shapiro_p_resid = tryCatch(shapiro.test(residuals(f))$p.value, error = function(e) NA))
}
tab <- table(ss$patient, ss$group)
pairs <- rownames(tab)[tab[, "normal"] == 1 & tab[, "tumour"] == 1]; stopifnot(length(pairs) == 8)   # 8 tumour/normal pairs (results/09_verification/sample_map_with_patient.tsv)
paired <- function(y) {
  d <- sapply(pairs, function(p) y[ss$patient == p & ss$group == "tumour"] - y[ss$patient == p & ss$group == "normal"])
  t <- t.test(d); w <- suppressWarnings(wilcox.test(d, exact = TRUE))
  c(n_pairs = length(d), mean_diff = mean(d), ci_low = t$conf.int[1], ci_high = t$conf.int[2], t_p = t$p.value, wilcoxon_p = w$p.value,
    n_pairs_same_sign_as_mean = sum(sign(d) == sign(mean(d))))
}
analyse <- function(ids, mat_re) {
  r <- t(sapply(ids, function(g) fit_gene(mat_re[g, ]))); r <- as.data.frame(r)
  r$p_BH <- p.adjust(r$p, "BH"); r$gene_id <- ids; r
}

dir.create(out, recursive = TRUE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- 0. reference-gene behaviour (circular: selected on these samples for exactly this property)
rs <- as.data.frame(t(fit_gene(refmean))); rs$what <- "mean log2 expression of the 10 reference genes"
rs$sd_across_samples <- sd(refmean)
w(rs[, c("what", "est", "ci_low", "ci_high", "p", "sd_across_samples")], "reference_mean_tumour_vs_normal.tsv")
ref_cor <- cor(t(log2(nc[ref_rows, ] + 1)))
w(data.frame(symbol = ref_sym, gene_id = ref_rows, mean_log2 = rowMeans(log2(nc[ref_rows, ] + 1)), sd_log2 = apply(log2(nc[ref_rows, ] + 1), 1, sd),
             mean_corr_with_other_refs = (rowSums(ref_cor) - 1) / 9), "reference_genes_summary.tsv")

# ---- 1. per-sample relative expression of the 9 model genes
mg <- model$gene_id
long <- do.call(rbind, lapply(mg, function(g) data.frame(sample_id = colnames(nc), patient = ss$patient, study = ss$study, group = ss$group,
  symbol = model$symbol[model$gene_id == g], gene_id = g, model_set = model$set[model$gene_id == g],
  log2_norm_count_pc1 = log2(nc[g, ] + 1), ref_mean_log2 = refmean, relative_expression = RE[g, ])))
w(long, "model9_relative_expression_per_sample.tsv")

# ---- 2. tumour vs normal, 9 model genes (primary + sensitivities)
m9 <- analyse(mg, RE)
m9$symbol <- model$symbol[match(m9$gene_id, model$gene_id)]; m9$model_set <- model$set[match(m9$gene_id, model$gene_id)]
pr <- as.data.frame(t(sapply(mg, function(g) paired(RE[g, ])))); names(pr) <- paste0("paired_", names(pr))
sens <- sapply(list(pc0.5 = re_of(nc, 0.5), pc4 = re_of(nc, 4), raw_counts_refs_only_pc1 = re_of(raw, 1)), function(M) sapply(mg, function(g) fit_gene(M[g, ])["est"]))
colnames(sens) <- paste0("est_", colnames(sens))
d <- de[match(mg, de$gene_id), ]
m9 <- cbind(m9[, c("symbol", "gene_id", "model_set", "est", "se", "ci_low", "ci_high", "p", "p_BH", "shapiro_p_resid")], pr, sens,
            deseq2_log2FC_unshrunk = d$log2FoldChange, deseq2_log2FC_shrunk = d$log2FoldChange_shrunk, deseq2_padj = d$padj)
m9$direction_matches_model_set <- ifelse(m9$model_set == "up", m9$est > 0, m9$est < 0)
w(m9, "model9_tumour_vs_normal.tsv")

# ---- 3. exploratory: the 1,390 candidates (primary model only)
cx <- analyse(cand$gene_id, RE)
cx$symbol <- cand$symbol; cx$direction_fc <- cand$direction_fc
cx$deseq2_log2FC_unshrunk <- de$log2FoldChange[match(cx$gene_id, de$gene_id)]
cx$same_sign_as_deseq2 <- sign(cx$est) == sign(cx$deseq2_log2FC_unshrunk)
w(cx[, c("symbol", "gene_id", "direction_fc", "est", "ci_low", "ci_high", "p", "p_BH", "deseq2_log2FC_unshrunk", "same_sign_as_deseq2")], "candidates1390_tumour_vs_normal_EXPLORATORY.tsv")

# ---- 4. summary for the log
cat(sprintf("reference mean: tumour - normal %.3f (95%% CI %.3f to %.3f), p %.3g\n", rs$est, rs$ci_low, rs$ci_high, rs$p))
cat(sprintf("9 model genes: %d of 9 with BH p < 0.05; direction matches model set for %d of 9; paired pairs n = %d\n",
            sum(m9$p_BH < 0.05), sum(m9$direction_matches_model_set), length(pairs)))
cat(sprintf("1,390 candidates: %d with BH p < 0.05; same sign as DESeq2 for %d; max |RE effect - DESeq2 unshrunk LFC| %.3f\n",
            sum(cx$p_BH < 0.05), sum(cx$same_sign_as_deseq2), max(abs(cx$est - cx$deseq2_log2FC_unshrunk))))
writeLines(c(paste("R", getRversion()), paste("stats (base R)"), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("command:", paste(commandArgs(), collapse = " "))), file.path(out, "tool_versions.txt"))
