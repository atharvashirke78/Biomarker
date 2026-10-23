#!/usr/bin/env Rscript
# PREREG (g): time ONE robustness refit (same DESeq2 + apeglm call as candidate_rules.R fit(), design ~ group, pre-filter >= 10 counts in >= n,
# n = smallest TRAINING group) on one 80 % patient-level subsample of the training set (seed 42). Writes elapsed seconds only; no gene results.
# usage: time_probe.R <counts.tsv> <clinical.tsv> <holdout_samples.txt> <out.txt>
suppressPackageStartupMessages({ library(DESeq2); library(apeglm) })
a <- commandArgs(trailingOnly = TRUE); if (file.exists(a[4])) stop(a[4], " exists")
cn <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
cl <- read.delim(a[2], stringsAsFactors = FALSE, colClasses = "character"); rownames(cl) <- cl$sample_id
ho <- readLines(a[3]); tr <- setdiff(colnames(cn), ho); ss <- cl[tr, ]; ss$group <- relevel(factor(ss$group), ref = "normal")
MIN_SAMPLES <- min(table(ss$group))
set.seed(42); pats <- unique(ss$patient_id); keep <- sample(pats, ceiling(0.80 * length(pats))); cols <- tr[ss$patient_id %in% keep]
t0 <- Sys.time()
d <- DESeqDataSetFromMatrix(round(cn[, cols]), colData = ss[cols, ], design = ~ group)
d <- d[rowSums(counts(d) >= 10) >= MIN_SAMPLES, ]; d <- DESeq(d, quiet = TRUE)
res <- results(d, name = "group_tumour_vs_normal", alpha = 0.05); shr <- lfcShrink(d, coef = "group_tumour_vs_normal", type = "apeglm", res = res, quiet = TRUE)
t <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
R <- if (20 * t <= 7200) 20 else max(5, floor(7200 / t))
writeLines(c(sprintf("probe samples %d (of %d training), genes after pre-filter %d, elapsed_seconds %.1f", length(cols), length(tr), nrow(d), t),
             sprintf("projected 20 refits = %.0f s (cap 7200 s) -> R = %d, pass >= %d%s", 20 * t, R, ceiling(0.90 * R - 1e-9), if (R * t > 7200) " (OVERRUN: even R = 5 exceeds 2 h)" else "")), a[4])
cat(readLines(a[4]), sep = "\n")
