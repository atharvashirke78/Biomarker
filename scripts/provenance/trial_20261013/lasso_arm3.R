#!/usr/bin/env Rscript
# LASSO trial, ARM 3 (reference-gene normalisation). Run from the project root with env/prefix/r_classifier/bin/Rscript, after lasso_stage1_deseq2.R (which supplies the per-fold features
# = rule 1 + baseMean >= 20 from the DESeq2 refit WITHOUT the test sample(s)). Same design and settings as lasso_stage2_glmnet.R (functions copied unchanged): featureCounts, seed 42 (set.seed(42 + fold index)),
# glmnet alpha = 1 binomial, lambda.1se (lambda.min also reported), inner 5-fold CV stratified by group, no study term, Youden cutoff on inner-CV out-of-fold scores, positive = tumour.
# Arm 3 value of gene g in sample s:  x = log2(count_gs + 1) - mean over the 10 reference genes r of log2(count_rs + 1)      (raw featureCounts counts; the subtraction removes library size; no training fit, no size factors)
# Reference genes: config/reference_genes.yaml (status PROPOSED), chosen using ALL 19 samples, so stability is partly built in. The 15 validation samples are not used.
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
s1 <- "results/14_lasso_trial/stage1"; out <- "results/14_lasso_trial/arm3"
sm <- read.delim("results/sample_map.tsv", stringsAsFactors = FALSE); rownames(sm) <- sm$srr; y_all <- as.integer(sm$group == "tumour"); names(y_all) <- sm$srr
sym <- read.delim("results/06_dge/deseq2_results_all.tsv", stringsAsFactors = FALSE); sym <- setNames(sym$gene_name, sym$gene_id)
man <- read.delim(file.path(s1, "fold_manifest.tsv"), stringsAsFactors = FALSE)
cn <- as.matrix(read.delim("results/04_counts/counts_featurecounts.tsv", row.names = 1, check.names = FALSE)); stopifnot(ncol(cn) == 19)
ref_ids <- sub("^.*gene_id:\\s*(ENSG[0-9]+).*$", "\\1", grep("gene_id: ENSG", readLines("config/reference_genes.yaml"), value = TRUE)); stopifnot(length(ref_ids) == 10)
vid <- sub("\\..*$", "", rownames(cn)); ri <- match(ref_ids, vid); stopifnot(!anyNA(ri))
L <- log2(cn + 1); refmean <- colMeans(L[ri, , drop = FALSE]); A3 <- sweep(L, 2, refmean, "-")                     # Arm 3 matrix, all genes x 19 samples
writeLines(c("reference genes:", paste(ref_ids, sym[rownames(cn)[ri]]), "", "Arm 3 value = log2(count + 1) - mean(log2(count + 1) of the 10 reference genes), per sample", "reference set chosen using all 19 samples (stability partly built in)"), file.path(out, "arm3_definition.txt"))

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
  list(lambda = lam, selected = sel, coefs = cf[sel, 1], score = sc, prob = plogis(sc), thr = thr, pred = pred, note = note)
}
preds <- list(); sels <- list()
for (lam_type in c("lambda.1se", "lambda.min")) for (i in seq_len(nrow(man))) {
  fold <- man$fold[i]; tr <- strsplit(man$train[i], ";")[[1]]; te <- strsplit(man$test[i], ";")[[1]]
  f <- intersect(readLines(file.path(s1, "features", paste0(fold, ".txt"))), rownames(cn)); fi <- match(f, rownames(cn))   # feature IDs are versioned, as the matrix rownames
  Xtr <- t(A3[fi, tr, drop = FALSE]); Xte <- t(A3[fi, te, drop = FALSE]); colnames(Xtr) <- colnames(Xte) <- rownames(cn)[fi]
  set.seed(42 + i); r <- fit_fold(Xtr, y_all[tr], Xte, lam_type)
  preds[[length(preds) + 1]] <- data.frame(arm = "arm3_reference_genes", lambda_type = lam_type, fold = fold, sample = te, label = sm[te, "label"], study = sm[te, "study"], truth = sm[te, "group"], y = y_all[te],
    n_candidate_features = ncol(Xtr), lambda = r$lambda, n_selected = length(r$selected), score_link = r$score, prob_tumour = r$prob, youden_cutoff_link = r$thr, pred_tumour = r$pred, note = r$note, stringsAsFactors = FALSE)
  if (length(r$selected)) sels[[length(sels) + 1]] <- data.frame(arm = "arm3_reference_genes", lambda_type = lam_type, fold = fold, gene_id = r$selected, symbol = sym[r$selected], coef = r$coefs, stringsAsFactors = FALSE)
}
P <- do.call(rbind, preds); S <- do.call(rbind, sels)
write.table(P, file.path(out, "fold_predictions_arm3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); write.table(S, file.path(out, "selected_genes_per_fold_arm3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summ <- function(d) { r <- roc(d$y, d$prob_tumour, levels = c(0, 1), direction = "<", quiet = TRUE); ci <- ci.auc(r, method = "delong")
  tp <- sum(d$pred_tumour == 1 & d$y == 1); fn <- sum(d$pred_tumour == 0 & d$y == 1); tn <- sum(d$pred_tumour == 0 & d$y == 0); fp <- sum(d$pred_tumour == 1 & d$y == 0)
  c(n_test = nrow(d), AUC = as.numeric(auc(r)), AUC_ci_low = ci[1], AUC_ci_high = ci[3], sensitivity = tp / (tp + fn), specificity = tn / (tn + fp), TP = tp, FN = fn, TN = tn, FP = fp) }
rows <- list()
for (lt in unique(P$lambda_type)) { d <- P[P$lambda_type == lt, ]; loo <- d[grepl("^loo_", d$fold), ]
  rows[[length(rows) + 1]] <- data.frame(design = "leave-one-out (pooled 19)", arm = "arm3_reference_genes", lambda_type = lt, t(summ(loo)), n_selected_median = median(loo$n_selected), n_selected_min = min(loo$n_selected), n_selected_max = max(loo$n_selected), n_candidate_features_median = median(loo$n_candidate_features), stringsAsFactors = FALSE)
  for (fo in c("los_A_to_B", "los_B_to_A")) { x <- d[d$fold == fo, ]; rows[[length(rows) + 1]] <- data.frame(design = paste("leave-one-study-out", fo), arm = "arm3_reference_genes", lambda_type = lt, t(summ(x)), n_selected_median = x$n_selected[1], n_selected_min = x$n_selected[1], n_selected_max = x$n_selected[1], n_candidate_features_median = x$n_candidate_features[1], stringsAsFactors = FALSE) } }
A3s <- do.call(rbind, rows); write.table(A3s, file.path(out, "auc_summary_arm3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
old <- read.delim("results/14_lasso_trial/auc_summary.tsv", stringsAsFactors = FALSE); cmp <- rbind(old, A3s[, names(old)]); write.table(cmp, file.path(out, "auc_summary_arms_1_2_3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(cmp[cmp$lambda_type == "lambda.1se", c("design", "arm", "AUC", "sensitivity", "specificity", "TP", "FN", "TN", "FP", "n_selected_median")], digits = 3)
for (lt in unique(S$lambda_type)) { s <- S[S$lambda_type == lt & grepl("^loo_", S$fold), ]; tb <- sort(table(s$gene_id), decreasing = TRUE)
  f <- data.frame(gene_id = names(tb), symbol = sym[names(tb)], n_folds_selected = as.integer(tb), fraction_of_19 = as.integer(tb) / 19, stringsAsFactors = FALSE)
  write.table(f, file.path(out, sprintf("selection_frequency_arm3_%s.tsv", lt)), sep = "\t", quote = FALSE, row.names = FALSE); cat("\narm3", lt, "top 10 most selected\n"); print(head(f, 10), row.names = FALSE) }

# ---- study shift of the genes selected in the leave-one-study-out models (train study vs test study), Arm 2 versus Arm 3
rd <- function(f) read.delim(f, row.names = 1, check.names = FALSE)
S2 <- read.delim("results/14_lasso_trial/selected_genes_per_fold.tsv", stringsAsFactors = FALSE)
sh <- list()
for (fo in c("los_A_to_B", "los_B_to_A")) {
  i <- which(man$fold == fo); tr <- strsplit(man$train[i], ";")[[1]]; te <- strsplit(man$test[i], ";")[[1]]
  M2tr <- rd(file.path(s1, "arm2", paste0(fo, "_train.tsv"))); M2te <- rd(file.path(s1, "arm2", paste0(fo, "_test.tsv")))
  sets <- list(arm2_selected = S2$gene_id[S2$arm == "arm2_train_only" & S2$lambda_type == "lambda.1se" & S2$fold == fo], arm3_selected = S$gene_id[S$lambda_type == "lambda.1se" & S$fold == fo])
  for (sn in names(sets)) { g <- sets[[sn]]
    for (arm in c("arm2_train_only", "arm3_reference_genes")) {
      if (arm == "arm2_train_only") { gg <- intersect(g, rownames(M2tr)); Xtr <- as.matrix(M2tr[gg, tr, drop = FALSE]); Xte <- as.matrix(M2te[gg, te, drop = FALSE]) }
      else { gg <- intersect(g, rownames(cn)); Xtr <- A3[gg, tr, drop = FALSE]; Xte <- A3[gg, te, drop = FALSE] }
      ntr <- tr[sm[tr, "group"] == "normal"]; nte <- te[sm[te, "group"] == "normal"]
      Xn_tr <- if (arm == "arm2_train_only") as.matrix(M2tr[gg, ntr, drop = FALSE]) else A3[gg, ntr, drop = FALSE]
      Xn_te <- if (arm == "arm2_train_only") as.matrix(M2te[gg, nte, drop = FALSE]) else A3[gg, nte, drop = FALSE]
      d_all <- apply(Xte, 1, median) - apply(Xtr, 1, median); d_nor <- apply(Xn_te, 1, median) - apply(Xn_tr, 1, median); sdtr <- apply(Xtr, 1, sd)
      sh[[length(sh) + 1]] <- data.frame(fold = fo, gene_set = sn, arm = arm, n_genes = length(gg), median_shift_all = median(d_all), median_abs_shift_all = median(abs(d_all)),
        median_abs_shift_normals_only = median(abs(d_nor)), median_abs_shift_in_train_SD = median(abs(d_all) / sdtr), stringsAsFactors = FALSE) } } }
SH <- do.call(rbind, sh); write.table(SH, file.path(out, "study_shift_selected_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(SH, digits = 3)
writeLines(c(paste("R", getRversion()), paste("glmnet", packageVersion("glmnet")), paste("pROC", packageVersion("pROC")), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions_arm3.txt"))
