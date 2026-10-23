#!/usr/bin/env Rscript
# Diagnostic model fit for -entry DISCOVER (Phase 1, 2026-10-15). COPY of the model part of results_rc1prep_19/lock_proposal_v2/fit_final_model_v2.R
# (the original is unchanged): same matrix (log2 CPM+1 of the featureCounts counts, CPM total = column sum), same features (candidates passing all
# applied rules, in candidate-table order), same glmnet call (binomial, alpha 1, foldid K folds stratified by group x batch, seed lasso_seed,
# type.measure deviance, lambda.1se, standardize = glmnet default), same cutoff (median of K inner-fold Youden cutoffs, seed cutoff_seed).
# Generalised: groups, baseline and batch from the DISCOVER sample sheet. Writes ONLY into the current folder (a new run folder), never into lock/.
# Not copied: the cross-study checks and the arm-4 selection-frequency column (trial artefacts of the lock proposal).
# usage: fit_diagnostic_model.R --counts <counts_training.tsv> --sheet <sample_sheet.discover.tsv> --candidates <candidate_table.tsv>
#        --de <deseq2_results_all.tsv> --questions config/questions.yaml --method-note <candidate_rules_method_note.txt>
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
a <- commandArgs(trailingOnly = TRUE); opt <- list()
for (i in seq(1, length(a), by = 2)) opt[[sub("^--", "", a[i])]] <- a[i + 1]
for (k in c("counts", "sheet", "candidates", "de", "questions", "method-note")) if (is.null(opt[[k]])) stop("missing --", k)
read_key <- function(file, key, sep) { l <- grep(paste0("^\\s*", key, "\\s*", sep), readLines(file), value = TRUE); stopifnot(length(l) == 1)
  trimws(gsub("['\"]", "", sub("(#|//).*$", "", sub(paste0("^\\s*", key, "\\s*", sep), "", l)))) }
SEED <- as.integer(read_key(opt$questions, "lasso_seed", ":")); CSEED <- as.integer(read_key(opt$questions, "cutoff_seed", ":")); K <- as.integer(read_key(opt$questions, "cv_folds", ":"))
cn <- as.matrix(read.delim(opt$counts, row.names = 1, check.names = FALSE))
ss <- read.delim(opt$sheet, stringsAsFactors = FALSE, colClasses = "character"); rownames(ss) <- ss$sample_id; ss <- ss[colnames(cn), ]
base_group <- unique(ss$group[toupper(ss$baseline) == "TRUE"]); stopifnot(length(base_group) == 1); pos_group <- setdiff(unique(ss$group), base_group)
y <- as.integer(ss$group == pos_group); names(y) <- colnames(cn); batch <- setNames(ss$batch, colnames(cn))
tot <- colSums(cn); L4 <- log2(sweep(cn, 2, tot, "/") * 1e6 + 1)                 # log2(CPM + 1), CPM total = featureCounts assigned total
res <- read.delim(opt$de, stringsAsFactors = FALSE)
ct <- read.delim(opt$candidates, stringsAsFactors = FALSE)
pool <- ct$gene_id[ct$all_rules_pass %in% c(TRUE, "TRUE")]                       # candidate-table order (= final_candidates order in the 19-sample run)
feat <- intersect(pool, rownames(cn)); cat("candidate features:", length(feat), "\n"); stopifnot(length(feat) >= 2)
strat_foldid <- function(strata, k = K) { id <- integer(length(strata)); for (cl in unique(strata)) { ix <- which(strata == cl); id[ix] <- sample(rep_len(seq_len(k), length(ix))) }; id }
youden_thr <- function(yy, s) { r <- roc(yy, s, levels = c(0, 1), direction = "<", quiet = TRUE); th <- coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[, 1]; median(th[is.finite(th)]) }
sc_fun <- function(up, dn, smp) sapply(smp, function(s) { u <- if (length(up)) mean(L4[up, s]) else 0; d <- if (length(dn)) mean(L4[dn, s]) else 0; u - d })
fit_model <- function(features, tr, strata_tr) {
  X <- t(L4[features, tr, drop = FALSE]); yt <- y[tr]
  set.seed(SEED); fid <- strat_foldid(strata_tr)
  cv <- cv.glmnet(X, yt, family = "binomial", alpha = 1, foldid = fid, type.measure = "deviance")
  cf <- as.matrix(coef(cv, s = "lambda.1se")); cf <- cf[rownames(cf) != "(Intercept)", 1]; cf <- cf[cf != 0]
  cfm <- as.matrix(coef(cv, s = "lambda.min")); cfm <- cfm[rownames(cfm) != "(Intercept)", 1]
  up <- names(cf)[cf > 0]; dn <- names(cf)[cf < 0]
  set.seed(CSEED); fid2 <- strat_foldid(strata_tr); thr_f <- numeric(K); s_tr <- sc_fun(up, dn, tr)
  for (k in seq_len(K)) { itr <- fid2 != k; thr_f[k] <- youden_thr(yt[itr], s_tr[itr]) }
  list(cf = cf, n_lambda_min = sum(cfm != 0), lambda = cv$lambda.1se, up = up, dn = dn, thr = median(thr_f), thr_folds = thr_f, thr_insample = youden_thr(yt, s_tr))
}
perf <- function(yy, s, thr) { r <- roc(yy, s, levels = c(0, 1), direction = "<", quiet = TRUE); ci <- suppressWarnings(ci.auc(r, method = "delong")); p <- as.integer(s > thr)
  data.frame(n = length(yy), AUC = as.numeric(auc(r)), ci_low = ci[1], ci_high = ci[3], TP = sum(p == 1 & yy == 1), FN = sum(p == 0 & yy == 1), TN = sum(p == 0 & yy == 0), FP = sum(p == 1 & yy == 0)) }
allS <- colnames(cn); strata <- paste(y[allS], batch[allS])
fm <- fit_model(feat, allS, strata)
if (length(fm$up) == 0 || length(fm$dn) == 0) message("NOTE: the selected genes are all in one direction (n_up ", length(fm$up), ", n_down ", length(fm$dn), "); the score uses 0 for the empty side")
cat("final: lambda.1se", fm$lambda, "| genes", length(fm$cf), "(up", length(fm$up), "down", length(fm$dn), ") | genes at lambda.min", fm$n_lambda_min, "\n")
cat("cutoff inner-fold:", signif(fm$thr_folds, 5), " median (final) =", fm$thr, " | in-sample Youden =", fm$thr_insample, "\n")
sAll <- sc_fun(fm$up, fm$dn, allS); P <- perf(y, sAll, fm$thr)
mod <- data.frame(gene_id = names(fm$cf), gene_id_unversioned = sub("\\..*$", "", names(fm$cf)), symbol = res$gene_name[match(names(fm$cf), res$gene_id)],
                  biotype = res$gene_type[match(names(fm$cf), res$gene_id)], set = ifelse(fm$cf > 0, "up", "down"), sign = ifelse(fm$cf > 0, "+1", "-1"),
                  lasso_coef_lambda1se = signif(fm$cf, 6), baseMean_fc = res$baseMean[match(names(fm$cf), res$gene_id)],
                  log2FC_shrunk_fc = res$log2FoldChange_shrunk[match(names(fm$cf), res$gene_id)], stringsAsFactors = FALSE)
mod <- mod[order(mod$set, mod$gene_id), ]
write.table(mod, "model_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
spec <- c("# DISCOVER diagnostic model (development fit; NOT the locked model lock/diagnostic_v1, which is never changed by DISCOVER)",
 "model: balanced_updown_score_log2CPM",
 "score(sample) = mean(log2(CPM + 1) over the UP genes) - mean(log2(CPM + 1) over the DOWN genes); CPM = featureCounts count / sample's featureCounts assigned total * 1e6",
 sprintf("call %s if score > cutoff (strictly greater); %s otherwise", pos_group, base_group),
 sprintf("n_up = %d, n_down = %d (genes and signs in model_genes.tsv)", length(fm$up), length(fm$dn)),
 sprintf("cutoff = %.10g   (median of the %d inner-fold Youden cutoffs: %s)", fm$thr, K, paste(signif(fm$thr_folds, 6), collapse = ", ")),
 sprintf("seed = %d (LASSO CV folds) and %d (inner cutoff folds); K = %d; folds stratified by group x batch", SEED, CSEED, K),
 sprintf("lasso: glmnet binomial alpha = 1, lambda.1se = %.8g (cv.glmnet deviance), candidate features = %d genes passing all applied candidate rules", fm$lambda, length(feat)),
 sprintf("training samples: %d (%s)", length(allS), paste(sprintf("%s %d", names(table(ss$group)), table(ss$group)), collapse = ", ")),
 sprintf("reference for comparison only (not the cutoff): in-sample Youden cutoff = %.10g", fm$thr_insample),
 paste("candidate rules:", readLines(opt[["method-note"]])))
writeLines(spec, "model_spec.txt")
write.table(data.frame(fold = seq_len(K), inner_fold_cutoff = fm$thr_folds), "inner_fold_cutoffs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = allS, batch = batch[allS], truth = ss$group, score = sAll, pred_positive = as.integer(sAll > fm$thr)), "development_scores.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
P$label <- "development data (training samples): OPTIMISTIC, genes and cutoff chosen on these samples; NOT validation"
write.table(P, "development_performance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("R", getRversion()), paste("glmnet", packageVersion("glmnet")), paste("pROC", packageVersion("pROC"))), "tool_versions.txt")
