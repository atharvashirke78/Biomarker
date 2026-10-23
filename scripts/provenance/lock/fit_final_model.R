#!/usr/bin/env Rscript
# DRAFT (step 12a). Final primary model = Arm 4 balanced score, fitted on the 19 DEVELOPMENT samples only (counts_featurecounts.tsv has 19 columns; the 15 validation samples are not in it).
# Run from project root: env/prefix/r_classifier/bin/Rscript results_rc1prep_19/lock_proposal/fit_final_model.R   (writes only into results_rc1prep_19/lock_proposal/)
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
SEED <- 42; K <- 5; out <- "results_rc1prep_19/lock_proposal"; R <- "results_rc1prep_19"
cn <- as.matrix(read.delim("results/04_counts/counts_featurecounts.tsv", row.names = 1, check.names = FALSE)); stopifnot(ncol(cn) == 19)
sm <- read.delim("results/sample_map.tsv", stringsAsFactors = FALSE); rownames(sm) <- sm$srr; stopifnot(all(colnames(cn) %in% sm$srr)); ss <- sm[colnames(cn), ]
y <- as.integer(ss$group == "tumour"); names(y) <- colnames(cn); study <- setNames(ss$study, colnames(cn))
tot <- colSums(cn); L4 <- log2(sweep(cn, 2, tot, "/") * 1e6 + 1)                 # log2(CPM + 1), CPM total = featureCounts assigned total (checked in the trial)
res <- read.delim(file.path(R, "deseq2_results_all.tsv"), stringsAsFactors = FALSE)   # NEW 19-sample fit, featureCounts route
feat19 <- res$gene_id[!is.na(res$padj) & res$padj < 0.05 & !is.na(res$log2FoldChange_shrunk) & abs(res$log2FoldChange_shrunk) >= 1 & res$baseMean >= 20]
feat19 <- intersect(feat19, rownames(cn)); cat("candidate features for the all-19 LASSO (rule 1 + baseMean >= 20, featureCounts):", length(feat19), "\n")
strat_foldid <- function(strata, k = K) { id <- integer(length(strata)); for (cl in unique(strata)) { ix <- which(strata == cl); id[ix] <- sample(rep_len(seq_len(k), length(ix))) }; id }
youden_thr <- function(yy, s) { r <- roc(yy, s, levels = c(0, 1), direction = "<", quiet = TRUE); th <- coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[, 1]; median(th[is.finite(th)]) }
sc_fun <- function(up, dn, smp) sapply(smp, function(s) { u <- if (length(up)) mean(L4[up, s]) else 0; d <- if (length(dn)) mean(L4[dn, s]) else 0; u - d })
fit_model <- function(features, tr, strata_tr, seed = SEED) {
  X <- t(L4[features, tr, drop = FALSE]); yt <- y[tr]
  set.seed(seed); fid <- strat_foldid(strata_tr)
  cv <- cv.glmnet(X, yt, family = "binomial", alpha = 1, foldid = fid, type.measure = "deviance")
  cf <- as.matrix(coef(cv, s = "lambda.1se")); cf <- cf[rownames(cf) != "(Intercept)", 1]; cf <- cf[cf != 0]
  cfm <- as.matrix(coef(cv, s = "lambda.min")); cfm <- cfm[rownames(cfm) != "(Intercept)", 1]
  up <- names(cf)[cf > 0]; dn <- names(cf)[cf < 0]
  # cutoff: inner CV (K folds stratified by strata), Youden J cutoff on the inner-training part of each fold, median of the K cutoffs
  set.seed(seed + 1); fid2 <- strat_foldid(strata_tr); thr_f <- numeric(K); s_tr <- sc_fun(up, dn, tr)
  for (k in seq_len(K)) { itr <- fid2 != k; thr_f[k] <- youden_thr(yt[itr], s_tr[itr]) }
  list(cf = cf, n_lambda_min = sum(cfm != 0), lambda = cv$lambda.1se, up = up, dn = dn, thr = median(thr_f), thr_folds = thr_f, thr_insample = youden_thr(yt, s_tr), s_tr = s_tr)
}
perf <- function(yy, s, thr) { r <- roc(yy, s, levels = c(0, 1), direction = "<", quiet = TRUE); ci <- suppressWarnings(ci.auc(r, method = "delong")); p <- as.integer(s > thr)
  data.frame(n = length(yy), AUC = as.numeric(auc(r)), ci_low = ci[1], ci_high = ci[3], TP = sum(p == 1 & yy == 1), FN = sum(p == 0 & yy == 1), TN = sum(p == 0 & yy == 0), FP = sum(p == 1 & yy == 0)) }
# ---- final model on all 19
all19 <- colnames(cn); strata <- paste(y[all19], study[all19])
fm <- fit_model(feat19, all19, strata)
cat("final: lambda.1se", fm$lambda, "| genes", length(fm$cf), "(up", length(fm$up), "down", length(fm$dn), ") | genes at lambda.min", fm$n_lambda_min, "\n")
cat("cutoff inner-fold:", signif(fm$thr_folds, 5), " median (final) =", fm$thr, " | in-sample Youden (arm 4 trial method) =", fm$thr_insample, "\n")
s19 <- sc_fun(fm$up, fm$dn, all19); P19 <- perf(y, s19, fm$thr)
sym <- setNames(res$gene_name, res$gene_id); typ <- setNames(res$gene_type, res$gene_id)
sel <- read.delim(file.path(R, "lasso_arm4/selection_frequency_arm4.tsv"), stringsAsFactors = FALSE)
mod <- data.frame(gene_id = names(fm$cf), gene_id_unversioned = sub("\\..*$", "", names(fm$cf)), symbol = sym[names(fm$cf)], biotype = typ[names(fm$cf)], set = ifelse(fm$cf > 0, "up", "down"), sign = ifelse(fm$cf > 0, "+1", "-1"), lasso_coef_lambda1se = signif(fm$cf, 6), stringsAsFactors = FALSE)
mod$loo_selection_freq_of_19 <- ifelse(mod$gene_id %in% sel$gene_id, sel$n_folds_selected[match(mod$gene_id, sel$gene_id)], 0)
cand <- read.delim(file.path(R, "final_candidates/final_candidates_18of19.tsv"), stringsAsFactors = FALSE); mod$in_final_candidates_1390 <- mod$gene_id %in% cand$gene_id
mod$baseMean_fc <- res$baseMean[match(mod$gene_id, res$gene_id)]; mod$log2FC_shrunk_fc <- res$log2FoldChange_shrunk[match(mod$gene_id, res$gene_id)]
mod <- mod[order(mod$set, -mod$loo_selection_freq_of_19), ]
write.table(mod, file.path(out, "final_model_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(mod[, c("symbol", "set", "lasso_coef_lambda1se", "loo_selection_freq_of_19", "in_final_candidates_1390")], row.names = FALSE)
spec <- c("# DRAFT final primary model (not locked). Balanced up/down score, Arm 4 design.",
 "model: balanced_updown_score_log2CPM",
 "score(sample) = mean(log2(CPM + 1) over the UP genes) - mean(log2(CPM + 1) over the DOWN genes); CPM = featureCounts count / sample's featureCounts assigned total * 1e6",
 "call tumour if score > cutoff (strictly greater); normal otherwise",
 sprintf("n_up = %d, n_down = %d (genes and signs in final_model_genes.tsv)", length(fm$up), length(fm$dn)),
 sprintf("cutoff = %.10g   (median of the %d inner-fold Youden cutoffs: %s)", fm$thr, K, paste(signif(fm$thr_folds, 6), collapse = ", ")),
 sprintf("seed = %d (LASSO CV folds) and %d (inner cutoff folds); K = %d; folds stratified by group x study", SEED, SEED + 1, K),
 sprintf("lasso: glmnet binomial alpha = 1, lambda.1se = %.8g (cv.glmnet deviance), candidate features = rule 1 + baseMean >= 20 on the featureCounts 19-sample fit (%d genes)", fm$lambda, length(feat19)),
 "training samples: the 19 development samples (results/04_counts/counts_featurecounts.tsv); the 15 validation samples were not used",
 sprintf("reference for comparison only (not the lock cutoff): in-sample Youden cutoff = %.10g", fm$thr_insample))
writeLines(spec, file.path(out, "final_model_spec.txt"))
# ---- development performance (OPTIMISTIC: genes chosen on the same 19 samples)
cat("\ndevelopment-data performance at final cutoff (optimistic):\n"); print(P19)
# ---- cross-study checks (same procedure, features from stage 1 'los' folds of the NEW fit)
xs <- list(); XS <- NULL
for (fo in c("los_A_to_B", "los_B_to_A")) {
  man <- read.delim(file.path(R, "lasso_stage1/fold_manifest.tsv"), stringsAsFactors = FALSE); m <- man[man$fold == fo, ]
  tr <- strsplit(m$train, ";")[[1]]; te <- strsplit(m$test, ";")[[1]]
  f <- intersect(readLines(file.path(R, "lasso_stage1/features", paste0(fo, ".txt"))), rownames(cn))
  fx <- fit_model(f, tr, as.character(y[tr]))
  ste <- sc_fun(fx$up, fx$dn, te); pe <- perf(y[te], ste, fx$thr)
  xs[[fo]] <- data.frame(fold = fo, n_train = length(tr), n_test = length(te), n_features = length(f), n_up = length(fx$up), n_down = length(fx$dn), cutoff = fx$thr, cutoff_insample_youden = fx$thr_insample, pe, stringsAsFactors = FALSE) }
XS <- do.call(rbind, xs); print(XS, digits = 3); write.table(XS, file.path(out, "cross_study_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = all19, study = study[all19], truth = ss$group, score = s19, pred_tumour = as.integer(s19 > fm$thr)), file.path(out, "development_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(P19, file.path(out, "development_performance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("R", getRversion()), paste("glmnet", packageVersion("glmnet")), paste("pROC", packageVersion("pROC")), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions_fit.txt"))
