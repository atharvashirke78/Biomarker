# rc1-prep copy (tracked) of results/14_lasso_trial/arm4/lasso_arm4.R; original left unchanged as the record of the old run. Args: no defaults; refuses a non-empty output folder unless --force.
#!/usr/bin/env Rscript
# LASSO trial, ARM 4 (balanced score) + normal-sample study shift for Arms 1-4. Run from the project root with env/prefix/r_classifier/bin/Rscript after lasso_stage1_deseq2.R.
# Same design as before: featureCounts, seed 42 (set.seed(42 + fold index)), per-fold features = rule 1 + baseMean >= 20 from the DESeq2 refit WITHOUT the test sample(s), glmnet alpha = 1 binomial,
# lambda.1se, inner 5-fold CV stratified by group, no study term, positive = tumour. The 15 validation samples are not used.
# Arm 4: LASSO input = log2(CPM + 1) of the fold's features (CPM = count / that sample's own featureCounts assigned total * 1e6, checked against the featureCounts .summary files).
#   Selected genes are split by coefficient sign into an up set (coef > 0) and a down set (coef < 0).  score(sample) = mean(log2(CPM+1) of up set) - mean(log2(CPM+1) of down set),
#   for every sample; nothing depends on the training fit except which genes are in the two sets. Cutoff = Youden on the TRAINING scores only (in-sample, because the genes were selected on
#   those samples); when several thresholds tie (perfect separation) the median of the tied thresholds is used. An empty set contributes 0 (score = mean(up) or -mean(down)); no genes -> training majority class.
#   Pooled leave-one-out AUC uses the raw scores of the 19 held-out samples.
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
args <- commandArgs(trailingOnly = TRUE); force <- "--force" %in% args; args <- setdiff(args, "--force")
usage <- "usage: %s %s [--force]   (no defaults)"
if (length(args) != 7) stop(sprintf(usage, "lasso_arm4.R", "<run_dir> <stage1_dir> <stage2_dir> <arm3_dir> <counts_dir> <ref_genes_yaml> <out_dir>"))
run_dir <- args[1]; s1 <- args[2]; stage2_dir <- args[3]; arm3_dir <- args[4]; counts_dir <- args[5]; ref_genes_yaml <- args[6]; out <- args[7]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) > 0 && !force) stop("ERROR: output folder ", out, " already contains files; nothing written. Use a new folder (or --force).")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sm <- read.delim(file.path(dirname(normalizePath(counts_dir)), "sample_map.tsv"), stringsAsFactors = FALSE); rownames(sm) <- sm$srr; y_all <- as.integer(sm$group == "tumour"); names(y_all) <- sm$srr
sym <- read.delim(file.path(run_dir, "deseq2_results_all.tsv"), stringsAsFactors = FALSE); sym <- setNames(sym$gene_name, sym$gene_id)
man <- read.delim(file.path(s1, "fold_manifest.tsv"), stringsAsFactors = FALSE)
cn <- as.matrix(read.delim(file.path(counts_dir, "counts_featurecounts.tsv"), row.names = 1, check.names = FALSE)); stopifnot(ncol(cn) == 19)
tot <- colSums(cn)
asg <- sapply(colnames(cn), function(s) { d <- read.delim(file.path(counts_dir, s, paste0(s, ".featureCounts.txt.summary")), row.names = 1, check.names = FALSE); d["Assigned", 1] })
stopifnot(all(abs(asg - tot) < 1)); message("CPM totals = featureCounts Assigned for all 19 samples")
L4 <- log2(sweep(cn, 2, tot, "/") * 1e6 + 1)                                  # log2(CPM + 1), all genes x 19
ref_ids <- sub("^.*gene_id:\\s*(ENSG[0-9]+).*$", "\\1", grep("gene_id: ENSG", readLines(ref_genes_yaml), value = TRUE)); vid <- sub("\\..*$", "", rownames(cn))
A3 <- sweep(log2(cn + 1), 2, colMeans(log2(cn + 1)[match(ref_ids, vid), , drop = FALSE]), "-")        # Arm 3 matrix (as in arm3/lasso_arm3.R)
a1 <- read.delim(file.path(s1, "arm1_vst_all19.tsv"), row.names = 1, check.names = FALSE); rd <- function(f) read.delim(f, row.names = 1, check.names = FALSE)

strat_foldid <- function(y, k = 5) { id <- integer(length(y)); for (cl in unique(y)) { ix <- which(y == cl); id[ix] <- sample(rep_len(seq_len(k), length(ix))) }; id }
youden_thr <- function(y, s) { r <- roc(y, s, levels = c(0, 1), direction = "<", quiet = TRUE); th <- coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[, 1]; median(th[is.finite(th)]) }
# ---- Arm 4 fold
arm4_fold <- function(fold, i, tr, te) {
  f <- intersect(readLines(file.path(s1, "features", paste0(fold, ".txt"))), rownames(cn)); Xtr <- t(L4[f, tr, drop = FALSE]); ytr <- y_all[tr]
  set.seed(42 + i); foldid <- strat_foldid(ytr)
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, foldid = foldid, type.measure = "deviance")
  cf <- as.matrix(coef(cv, s = "lambda.1se")); cf <- cf[rownames(cf) != "(Intercept)", 1]; cf <- cf[cf != 0]
  up <- names(cf)[cf > 0]; dn <- names(cf)[cf < 0]
  sc <- function(s) { u <- if (length(up)) mean(L4[up, s]) else 0; d <- if (length(dn)) mean(L4[dn, s]) else 0; u - d }
  str <- sapply(tr, sc); ste <- sapply(te, sc)
  thr <- if (length(cf) == 0) NA_real_ else youden_thr(ytr, str)
  pred <- if (is.na(thr)) rep(as.integer(mean(ytr) >= 0.5), length(te)) else as.integer(ste > thr)
  list(up = up, dn = dn, cf = cf, str = str, ste = ste, thr = thr, pred = pred, n_cand = length(f))
}
P <- list(); S <- list(); TR <- list()
for (i in seq_len(nrow(man))) { fold <- man$fold[i]; tr <- strsplit(man$train[i], ";")[[1]]; te <- strsplit(man$test[i], ";")[[1]]
  r <- arm4_fold(fold, i, tr, te)
  P[[i]] <- data.frame(arm = "arm4_balanced_score", fold = fold, sample = te, label = sm[te, "label"], study = sm[te, "study"], truth = sm[te, "group"], y = y_all[te], n_candidate_features = r$n_cand,
    n_selected = length(r$cf), n_up = length(r$up), n_down = length(r$dn), score = r$ste, cutoff = r$thr, pred_tumour = r$pred, stringsAsFactors = FALSE)
  if (length(r$cf)) S[[i]] <- data.frame(arm = "arm4_balanced_score", fold = fold, gene_id = names(r$cf), symbol = sym[names(r$cf)], coef = r$cf, set = ifelse(r$cf > 0, "up", "down"), stringsAsFactors = FALSE)
  TR[[fold]] <- list(tr = tr, te = te, str = r$str, ste = r$ste) }
P <- do.call(rbind, P); S <- do.call(rbind, S)
write.table(P, file.path(out, "fold_predictions_arm4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); write.table(S, file.path(out, "selected_genes_per_fold_arm4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summ <- function(d, sc) { r <- roc(d$y, d[[sc]], levels = c(0, 1), direction = "<", quiet = TRUE); ci <- suppressWarnings(ci.auc(r, method = "delong"))
  tp <- sum(d$pred_tumour == 1 & d$y == 1); fn <- sum(d$pred_tumour == 0 & d$y == 1); tn <- sum(d$pred_tumour == 0 & d$y == 0); fp <- sum(d$pred_tumour == 1 & d$y == 0)
  data.frame(n_test = nrow(d), AUC = as.numeric(auc(r)), AUC_ci_low = ci[1], AUC_ci_high = ci[3], sensitivity = tp / (tp + fn), specificity = tn / (tn + fp), TP = tp, FN = fn, TN = tn, FP = fp) }
loo <- P[grepl("^loo_", P$fold), ]
A4 <- rbind(cbind(design = "leave-one-out (pooled 19)", arm = "arm4_balanced_score", summ(loo, "score"), n_selected_median = median(loo$n_selected), n_selected_min = min(loo$n_selected), n_selected_max = max(loo$n_selected)),
  do.call(rbind, lapply(c("los_A_to_B", "los_B_to_A"), function(fo) { x <- P[P$fold == fo, ]; cbind(design = paste("leave-one-study-out", fo), arm = "arm4_balanced_score", summ(x, "score"), n_selected_median = x$n_selected[1], n_selected_min = x$n_selected[1], n_selected_max = x$n_selected[1]) })))
write.table(A4, file.path(out, "auc_summary_arm4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(A4, digits = 3)
tb <- sort(table(S$gene_id[grepl("^loo_", S$fold)]), decreasing = TRUE); f4 <- data.frame(gene_id = names(tb), symbol = sym[names(tb)], n_folds_selected = as.integer(tb), fraction_of_19 = as.integer(tb) / 19)
write.table(f4, file.path(out, "selection_frequency_arm4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(head(f4, 10), row.names = FALSE)

# ---- normal-sample study shift, Arms 1-4, leave-one-study-out: refit Arms 1-3 (same seeds) to get link scores of training AND test samples (checked against the stored test scores)
link_scores <- function(arm, fold, i, tr, te) {
  feat <- readLines(file.path(s1, "features", paste0(fold, ".txt")))
  if (arm == "arm1_leaky") { f <- intersect(feat, rownames(a1)); Xtr <- t(a1[f, tr, drop = FALSE]); Xte <- t(a1[f, te, drop = FALSE]) }
  else if (arm == "arm2_train_only") { A <- rd(file.path(s1, "arm2", paste0(fold, "_train.tsv"))); B <- rd(file.path(s1, "arm2", paste0(fold, "_test.tsv"))); f <- intersect(rownames(A), rownames(B)); Xtr <- t(A[f, tr, drop = FALSE]); Xte <- t(B[f, te, drop = FALSE]) }
  else { f <- intersect(feat, rownames(cn)); Xtr <- t(A3[f, tr, drop = FALSE]); Xte <- t(A3[f, te, drop = FALSE]) }
  set.seed(42 + i); ytr <- y_all[tr]; cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, foldid = strat_foldid(ytr), type.measure = "deviance", keep = TRUE)
  list(tr = setNames(as.numeric(predict(cv, newx = Xtr, s = "lambda.1se", type = "link")), tr), te = setNames(as.numeric(predict(cv, newx = Xte, s = "lambda.1se", type = "link")), te)) }
stored <- rbind(read.delim(file.path(stage2_dir, "fold_predictions.tsv"), stringsAsFactors = FALSE), read.delim(file.path(arm3_dir, "fold_predictions_arm3.tsv"), stringsAsFactors = FALSE))
sh <- list()
shift_row <- function(arm, fo, tr, te, str, ste) {
  ntr <- tr[sm[tr, "group"] == "normal"]; nte <- te[sm[te, "group"] == "normal"]; ttr <- tr[sm[tr, "group"] == "tumour"]
  shift <- median(ste[nte]) - median(str[ntr]); gap <- median(str[ttr]) - median(str[ntr])
  data.frame(fold = fo, arm = arm, normal_median_train = median(str[ntr]), normal_median_test = median(ste[nte]), shift_normals = shift, train_class_gap = gap, shift_as_fraction_of_gap = shift / gap, stringsAsFactors = FALSE) }
for (fo in c("los_A_to_B", "los_B_to_A")) { i <- which(man$fold == fo); tr <- strsplit(man$train[i], ";")[[1]]; te <- strsplit(man$test[i], ";")[[1]]
  for (arm in c("arm1_leaky", "arm2_train_only", "arm3_reference_genes")) { ls <- link_scores(arm, fo, i, tr, te)
    ref <- stored[stored$arm == arm & stored$lambda_type == "lambda.1se" & stored$fold == fo, ]; dif <- max(abs(ref$score_link[match(te, ref$sample)] - ls$te)); message(arm, " ", fo, ": refit vs stored test score max abs diff ", signif(dif, 3)); stopifnot(dif < 1e-6)
    sh[[length(sh) + 1]] <- shift_row(arm, fo, tr, te, ls$tr, ls$te) }
  sh[[length(sh) + 1]] <- shift_row("arm4_balanced_score", fo, tr, te, TR[[fo]]$str, TR[[fo]]$ste) }
SH <- do.call(rbind, sh); write.table(SH, file.path(out, "normal_study_shift_all_arms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(SH, digits = 3)
writeLines(c(paste("R", getRversion()), paste("glmnet", packageVersion("glmnet")), paste("pROC", packageVersion("pROC")), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions_arm4.txt"))
