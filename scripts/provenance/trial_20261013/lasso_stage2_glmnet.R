#!/usr/bin/env Rscript
# LASSO trial, stage 2 (env/prefix/r_classifier: glmnet + pROC). Run from the project root after lasso_stage1_deseq2.R.
# For each fold: glmnet alpha = 1, binomial, on the TRAINING samples only; inner 5-fold CV stratified by group (seed 42), lambda.1se (primary) and lambda.min.
# No study term in the classifier. Cutoff = Youden on the INNER-CV out-of-fold scores (cv.glmnet keep = TRUE, fit.preval) at the chosen lambda; never the test sample.
# Positive class = tumour. Pooled leave-one-out ROC: held-out probabilities of the 19 folds (each from a different model) -> one pROC ROC, DeLong 95% CI.
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
set.seed(42)
s1 <- "results/14_lasso_trial/stage1"; out <- "results/14_lasso_trial"
sm <- read.delim("results/sample_map.tsv", stringsAsFactors = FALSE); rownames(sm) <- sm$srr; y_all <- as.integer(sm$group == "tumour"); names(y_all) <- sm$srr
sym <- read.delim("results/06_dge/deseq2_results_all.tsv", stringsAsFactors = FALSE); sym <- setNames(sym$gene_name, sym$gene_id)
man <- read.delim(file.path(s1, "fold_manifest.tsv"), stringsAsFactors = FALSE)
a1 <- read.delim(file.path(s1, "arm1_vst_all19.tsv"), row.names = 1, check.names = FALSE)
rd <- function(f) read.delim(f, row.names = 1, check.names = FALSE)

strat_foldid <- function(y, k = 5) { id <- integer(length(y)); for (cl in unique(y)) { ix <- which(y == cl); id[ix] <- sample(rep_len(seq_len(k), length(ix))) }; id }
fit_fold <- function(Xtr, ytr, Xte, lam_type) {
  foldid <- strat_foldid(ytr)
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, foldid = foldid, type.measure = "deviance", keep = TRUE)
  lam <- cv[[lam_type]]; li <- which(cv$lambda == lam)
  cf <- as.matrix(coef(cv, s = lam)); sel <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
  inner <- cv$fit.preval[, li]; note <- ""
  thr <- NA_real_
  if (length(sel) == 0) note <- "no genes selected: intercept-only, class = training majority"
  else if (anyNA(inner)) note <- "NA in inner-CV scores"
  else { r <- roc(ytr, inner, levels = c(0, 1), direction = "<", quiet = TRUE); thr <- as.numeric(coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1, 1]) }
  sc <- as.numeric(predict(cv, newx = Xte, s = lam, type = "link"))
  pred <- if (is.na(thr)) rep(as.integer(mean(ytr) >= 0.5), length(sc)) else as.integer(sc > thr)
  list(lambda = lam, selected = sel, coefs = cf[sel, 1], score = sc, prob = plogis(sc), thr = thr, pred = pred, note = note, n_inner_auc = if (!is.na(thr)) as.numeric(auc(roc(ytr, inner, levels = c(0, 1), direction = "<", quiet = TRUE))) else NA)
}
get_X <- function(arm, fold, tr, te) {
  feat <- readLines(file.path(s1, "features", paste0(fold, ".txt")))
  if (arm == "arm1_leaky") { f <- intersect(feat, rownames(a1)); list(Xtr = t(a1[f, tr, drop = FALSE]), Xte = t(a1[f, te, drop = FALSE])) }
  else { A <- rd(file.path(s1, "arm2", paste0(fold, "_train.tsv"))); B <- rd(file.path(s1, "arm2", paste0(fold, "_test.tsv"))); f <- intersect(rownames(A), rownames(B)); list(Xtr = t(A[f, tr, drop = FALSE]), Xte = t(B[f, te, drop = FALSE])) }
}
preds <- list(); sels <- list()
for (arm in c("arm1_leaky", "arm2_train_only")) for (lam_type in c("lambda.1se", "lambda.min")) for (i in seq_len(nrow(man))) {
  fold <- man$fold[i]; tr <- strsplit(man$train[i], ";")[[1]]; te <- strsplit(man$test[i], ";")[[1]]
  xx <- get_X(arm, fold, tr, te)
  set.seed(42 + i)                                                          # one reproducible stream per fold
  r <- fit_fold(xx$Xtr, y_all[tr], xx$Xte, lam_type)
  preds[[length(preds) + 1]] <- data.frame(arm = arm, lambda_type = lam_type, fold = fold, sample = te, label = sm[te, "label"], study = sm[te, "study"], truth = sm[te, "group"], y = y_all[te],
    n_candidate_features = ncol(xx$Xtr), lambda = r$lambda, n_selected = length(r$selected), score_link = r$score, prob_tumour = r$prob, youden_cutoff_link = r$thr, pred_tumour = r$pred, note = r$note, stringsAsFactors = FALSE)
  if (length(r$selected)) sels[[length(sels) + 1]] <- data.frame(arm = arm, lambda_type = lam_type, fold = fold, gene_id = r$selected, symbol = sym[r$selected], coef = r$coefs, stringsAsFactors = FALSE)
}
P <- do.call(rbind, preds); S <- do.call(rbind, sels)
write.table(P, file.path(out, "fold_predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); write.table(S, file.path(out, "selected_genes_per_fold.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summ <- function(d) { r <- roc(d$y, d$prob_tumour, levels = c(0, 1), direction = "<", quiet = TRUE); ci <- ci.auc(r, method = "delong")
  tp <- sum(d$pred_tumour == 1 & d$y == 1); fn <- sum(d$pred_tumour == 0 & d$y == 1); tn <- sum(d$pred_tumour == 0 & d$y == 0); fp <- sum(d$pred_tumour == 1 & d$y == 0)
  c(n_test = nrow(d), AUC = as.numeric(auc(r)), AUC_ci_low = ci[1], AUC_ci_high = ci[3], sensitivity = tp / (tp + fn), specificity = tn / (tn + fp), TP = tp, FN = fn, TN = tn, FP = fp) }
rows <- list()
for (arm in unique(P$arm)) for (lt in unique(P$lambda_type)) {
  d <- P[P$arm == arm & P$lambda_type == lt, ]
  loo <- d[grepl("^loo_", d$fold), ]; rows[[length(rows) + 1]] <- data.frame(design = "leave-one-out (pooled 19)", arm = arm, lambda_type = lt, t(summ(loo)),
    n_selected_median = median(loo$n_selected), n_selected_min = min(loo$n_selected), n_selected_max = max(loo$n_selected), n_candidate_features_median = median(loo$n_candidate_features), stringsAsFactors = FALSE)
  for (fo in c("los_A_to_B", "los_B_to_A")) { x <- d[d$fold == fo, ]; rows[[length(rows) + 1]] <- data.frame(design = paste("leave-one-study-out", fo), arm = arm, lambda_type = lt, t(summ(x)),
    n_selected_median = x$n_selected[1], n_selected_min = x$n_selected[1], n_selected_max = x$n_selected[1], n_candidate_features_median = x$n_candidate_features[1], stringsAsFactors = FALSE) }
}
A <- do.call(rbind, rows); write.table(A, file.path(out, "auc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(A, digits = 3)
for (arm in unique(S$arm)) for (lt in unique(S$lambda_type)) {
  s <- S[S$arm == arm & S$lambda_type == lt & grepl("^loo_", S$fold), ]; tb <- sort(table(s$gene_id), decreasing = TRUE)
  f <- data.frame(gene_id = names(tb), symbol = sym[names(tb)], n_folds_selected = as.integer(tb), fraction_of_19 = as.integer(tb) / 19, stringsAsFactors = FALSE)
  write.table(f, file.path(out, sprintf("selection_frequency_%s_%s.tsv", arm, lt)), sep = "\t", quote = FALSE, row.names = FALSE); cat("\n", arm, lt, "top 10 most selected (19 LOO folds)\n"); print(head(f, 10), row.names = FALSE)
}
writeLines(c(paste("R", getRversion()), paste("glmnet", packageVersion("glmnet")), paste("pROC", packageVersion("pROC")), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions_stage2.txt"))
