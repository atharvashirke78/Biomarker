#!/usr/bin/env Rscript
# Internal validation of the DISCOVER diagnostic model (upgrade H3c/H3d/H3e/H3f). ADDITIVE and REPORT-ONLY: it never feeds back into the model,
# never changes fit_diagnostic_model.R or any existing output. It re-implements the model-building recipe of fit_diagnostic_model.R as a function
# (log2(CPM+1); glmnet binomial alpha 1, foldid K folds stratified by group x batch, seed lasso_seed, deviance; lambda.1se; up/down sets from the
# coefficient signs; score = mean(up) - mean(down); cutoff = median of K inner-fold Youden thresholds, seed cutoff_seed) and repeats it ENTIRELY
# inside every resample: repeated nested CV, Harrell bootstrap optimism, bootstrap cutoff interval, per-gene selection frequency (subsamples),
# and the same recipe with lambda.min as a sensitivity row.
# LIMITATION (stated in every output): the candidate list (rules 1-3b) was selected on ALL samples and is NOT repeated inside the resamples.
# Only the LASSO + score + cutoff stage is validated; estimates remain optimistic. This is NOT independent validation.
# usage: internal_validation.R --counts <counts_training.tsv> --sheet <sample_sheet.discover.tsv> --candidates <candidate_table.tsv>
#        --genes <model_genes.tsv> --spec <model_spec.txt> --scores <development_scores.tsv> [--out .] [--repeats 20] [--bootstrap 200]
#        [--subsamples 100] [--check-only true]
# Always exits 0 (a failure becomes a 'NOT RUN: reason' row), so it can never fail the pipeline.
a <- commandArgs(trailingOnly = TRUE); opt <- list()
for (i in seq(1, length(a), by = 2)) opt[[sub("^--", "", a[i])]] <- a[i + 1]
OUT <- if (is.null(opt$out)) "." else opt$out
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NREP <- as.integer(if (is.null(opt$repeats)) 20 else opt$repeats)
NBOOT <- as.integer(if (is.null(opt$bootstrap)) 200 else opt$bootstrap)
NSUB <- as.integer(if (is.null(opt$subsamples)) 100 else opt$subsamples)
OUTER_SEED0 <- 1000L; BOOT_SEED0 <- 2000L; SUB_SEED0 <- 3000L        # fixed seeds: repeat r uses 1000 + r, bootstrap b 2000 + b, subsample s 3000 + s
LIMIT <- paste("LIMITATION: the candidate list (candidate rules 1-3b) was selected on ALL samples and the candidate screening stays OUTSIDE the resampling loop.",
               "Only the LASSO + score + cutoff stage is repeated inside each resample, so these estimates are still optimistic.",
               "This is internal validation of the existing model, NOT independent validation.")
f_sum <- file.path(OUT, "internal_validation_summary.tsv"); f_freq <- file.path(OUT, "selection_frequency.tsv"); f_md <- file.path(OUT, "INTERNAL_VALIDATION.md")
wtab <- function(d, f) write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
not_run <- function(reason) {
  lab <- paste0("NOT RUN: ", reason)
  wtab(data.frame(metric = "internal_validation", value = NA, ci_low = NA, ci_high = NA, label = lab), f_sum)
  wtab(data.frame(gene_id = lab, symbol = NA, direction = NA, frequency = NA, in_final_model = NA), f_freq)
  writeLines(c("# Internal validation of the diagnostic model", "", paste0("**", lab, "**"), "",
               "No internal-validation estimate was produced. The diagnostic model and every other pipeline output are unaffected.", "", LIMIT), f_md)
  cat(lab, "\n"); quit(save = "no", status = 0)
}
res_all <- tryCatch({
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
options(warn = -1)      # glmnet warns "class has fewer than 8 observations" in every small fit; the message is expected for small cohorts
for (k in c("counts", "sheet", "candidates", "genes", "spec", "scores")) if (is.null(opt[[k]])) stop("missing --", k)
spec <- readLines(opt$spec)
if (any(grepl("NOT FITTED|NO MODEL FITTED", spec[1:min(3, length(spec))]))) not_run("no model (model_spec.txt says NOT FITTED)")
gt <- read.delim(opt$genes, stringsAsFactors = FALSE)
if (nrow(gt) == 0) not_run("no model (model_genes.tsv is empty)")
# ---- settings of the pipeline's own fit, read from model_spec.txt (seeds, K, cutoff, patient-level folds)
sl <- grep("^seed = ", spec, value = TRUE); stopifnot(length(sl) == 1)
SEED <- as.integer(sub("^seed = (\\d+) \\(LASSO.*$", "\\1", sl)); CSEED <- as.integer(sub("^seed = \\d+ \\(LASSO CV folds\\) and (\\d+) .*$", "\\1", sl))
K <- as.integer(sub("^.*K = (\\d+);.*$", "\\1", sl)); PLC <- grepl("samples of one patient are kept together", sl)
SP_CUT <- as.numeric(sub("^cutoff = ([^ ]+) .*$", "\\1", grep("^cutoff = ", spec, value = TRUE)))
SP_LAM <- as.numeric(sub("^.*lambda.1se = ([^ ]+) \\(.*$", "\\1", grep("^lasso: ", spec, value = TRUE)))
cn <- as.matrix(read.delim(opt$counts, row.names = 1, check.names = FALSE))
ss <- read.delim(opt$sheet, stringsAsFactors = FALSE, colClasses = "character"); rownames(ss) <- ss$sample_id; ss <- ss[colnames(cn), ]
base_group <- unique(ss$group[toupper(ss$baseline) == "TRUE"]); stopifnot(length(base_group) == 1); pos_group <- setdiff(unique(ss$group), base_group)
y <- as.integer(ss$group == pos_group); names(y) <- colnames(cn); batch <- setNames(ss$batch, colnames(cn))
N <- ncol(cn)
if (N < 10) not_run(sprintf("only %d samples (at least 10 needed)", N))
if (min(table(y)) < 5) not_run(sprintf("a class has fewer than 5 samples (%s)", paste(sprintf("%s %d", c(base_group, pos_group), c(sum(y == 0), sum(y == 1))), collapse = ", ")))
pat <- setNames(if ("patient_id" %in% names(ss)) ss$patient_id else colnames(cn), colnames(cn))
L4 <- log2(sweep(cn, 2, colSums(cn), "/") * 1e6 + 1)
ct <- read.delim(opt$candidates, stringsAsFactors = FALSE)
pool <- ct$gene_id[ct$all_rules_pass %in% c(TRUE, "TRUE")]; feat <- intersect(pool, rownames(cn))
if (length(feat) < 2) not_run("fewer than 2 candidate genes")
sym <- setNames(ct$symbol, ct$gene_id)
REPEATED_PATIENTS <- PLC && any(duplicated(pat))
# ---- the recipe (same code and same random draws as fit_diagnostic_model.R; idx = column indices of L4, duplicates allowed in bootstrap samples)
strat_foldid <- function(strata, k = K) { id <- integer(length(strata)); for (cl in unique(strata)) { ix <- which(strata == cl); id[ix] <- sample(rep_len(seq_len(k), length(ix))) }; id }
strat_foldid_pat <- function(strata, patient, k = K) {
  up <- unique(patient); sig <- vapply(up, function(p) paste(sort(unique(strata[patient == p])), collapse = "|"), "")
  pid <- integer(length(up)); names(pid) <- up
  for (cl in unique(sig)) { ix <- which(sig == cl); pid[ix] <- sample(rep_len(seq_len(k), length(ix))) }
  unname(pid[as.character(patient)]) }
foldid <- function(strata_vec, patient_vec) if (REPEATED_PATIENTS) strat_foldid_pat(strata_vec, patient_vec) else strat_foldid(strata_vec)
youden_thr <- function(yy, s) { r <- roc(yy, s, levels = c(0, 1), direction = "<", quiet = TRUE); th <- coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[, 1]; median(th[is.finite(th)]) }
score_idx <- function(up, dn, idx) sapply(idx, function(j) { u <- if (length(up)) mean(L4[up, j]) else 0; d <- if (length(dn)) mean(L4[dn, j]) else 0; u - d })
cutoff_for <- function(up, dn, idx, yt, strata_tr, patt) {
  set.seed(CSEED); fid2 <- foldid(strata_tr, patt); thr_f <- numeric(K); s_tr <- score_idx(up, dn, idx)
  for (k in seq_len(K)) { itr <- fid2 != k; thr_f[k] <- youden_thr(yt[itr], s_tr[itr]) }
  list(thr = median(thr_f), thr_folds = thr_f) }
# fit_recipe: lambda.1se (the pipeline's model) and lambda.min (sensitivity), both with the full cutoff rule
fit_recipe <- function(idx, with_min = TRUE) {
  strata_tr <- paste(y[idx], batch[idx]); yt <- y[idx]; patt <- pat[idx]
  X <- t(L4[feat, idx, drop = FALSE])
  set.seed(SEED); fid <- foldid(strata_tr, patt)
  cv <- cv.glmnet(X, yt, family = "binomial", alpha = 1, foldid = fid, type.measure = "deviance")
  mk <- function(s) { cf <- as.matrix(coef(cv, s = s)); cf <- cf[rownames(cf) != "(Intercept)", 1]; cf <- cf[cf != 0]
    up <- names(cf)[cf > 0]; dn <- names(cf)[cf < 0]; ct_ <- cutoff_for(up, dn, idx, yt, strata_tr, patt)
    list(cf = cf, up = up, dn = dn, thr = ct_$thr, thr_folds = ct_$thr_folds) }
  r <- list(se = mk("lambda.1se"), lambda = cv$lambda.1se)
  if (with_min) r$mn <- mk("lambda.min")
  r }
auc_of <- function(yy, s) as.numeric(auc(roc(yy, s, levels = c(0, 1), direction = "<", quiet = TRUE)))
safe_fit <- function(idx) tryCatch(fit_recipe(idx), error = function(e) NULL)
sens_spec <- function(yy, s, thr) { p <- as.integer(s > thr); c(sens = sum(p == 1 & yy == 1) / sum(yy == 1), spec = sum(p == 0 & yy == 0) / sum(yy == 0)) }
ALL <- seq_len(N); strata_all <- paste(y, batch)
# ---- (1) reproduce the pipeline's model on all samples (all candidates in candidate-table order)
fm <- fit_recipe(ALL)
sAll <- score_idx(fm$se$up, fm$se$dn, ALL)
dv <- read.delim(opt$scores, stringsAsFactors = FALSE)
genes_ok <- setequal(names(fm$se$cf), gt$gene_id)
sc_ok <- all(dv$sample == colnames(cn)) && max(abs(dv$score - sAll)) < 1e-8
cut_ok <- abs(fm$se$thr - SP_CUT) < 1e-8; lam_ok <- abs(fm$lambda - SP_LAM) < 1e-7 * max(1, abs(SP_LAM)) || abs(fm$lambda - SP_LAM) < 5e-9
repro <- genes_ok && sc_ok && cut_ok && lam_ok
cat(sprintf("reproduction of the pipeline model: genes %s, scores %s, cutoff %s (%.10g vs %.10g), lambda.1se %s (%.15g vs %.10g in model_spec.txt)\n", genes_ok, sc_ok, cut_ok, fm$se$thr, SP_CUT, lam_ok, fm$lambda, SP_LAM))
if (!is.null(opt[["check-only"]])) { quit(save = "no", status = if (repro) 0 else 1) }
# ---- apparent performance of the pipeline's own model (development_scores.tsv)
dy <- as.integer(dv$truth == pos_group); app_auc <- auc_of(dy, dv$score)
app_se <- sum(dv$pred_positive == 1 & dy == 1) / sum(dy == 1); app_sp <- sum(dv$pred_positive == 0 & dy == 0) / sum(dy == 0)
# ---- (2) repeated nested CV: outer K-fold stratified by group x batch, R repeats; whole recipe on each outer training set
nested <- function(variant) {
  out <- data.frame(rep = seq_len(NREP), auc = NA_real_, sens = NA_real_, spec = NA_real_, failed = 0L)
  for (r in seq_len(NREP)) {
    set.seed(OUTER_SEED0 + r); ofid <- foldid(strata_all, pat); sc <- rep(NA_real_, N); pr <- rep(NA_integer_, N)
    for (k in seq_len(K)) {
      te <- which(ofid == k); tr <- which(ofid != k); m <- if (variant == "se") nested_fits[[r]][[k]]$se else nested_fits[[r]][[k]]$mn
      if (is.null(nested_fits[[r]][[k]])) { out$failed[r] <- out$failed[r] + 1L; next }
      s <- score_idx(m$up, m$dn, te); sc[te] <- s; pr[te] <- as.integer(s > m$thr) }
    ok <- !is.na(sc)
    if (sum(ok) > 1 && length(unique(y[ok])) == 2) { out$auc[r] <- auc_of(y[ok], sc[ok]); out$sens[r] <- sum(pr[ok] == 1 & y[ok] == 1) / sum(y[ok] == 1); out$spec[r] <- sum(pr[ok] == 0 & y[ok] == 0) / sum(y[ok] == 0) } }
  out }
nested_fits <- vector("list", NREP)
for (r in seq_len(NREP)) { set.seed(OUTER_SEED0 + r); ofid <- foldid(strata_all, pat); nested_fits[[r]] <- lapply(seq_len(K), function(k) safe_fit(which(ofid != k))) }
ncv_se <- nested("se"); ncv_mn <- nested("mn")
# ---- (3) Harrell bootstrap: resample within class (so both classes always present); whole recipe in each resample
boot <- vector("list", NBOOT)
for (b in seq_len(NBOOT)) {
  set.seed(BOOT_SEED0 + b); idx <- unlist(lapply(c(0, 1), function(cl) { ix <- which(y == cl); ix[sample.int(length(ix), length(ix), replace = TRUE)] }))
  m <- safe_fit(idx)
  boot[[b]] <- if (is.null(m)) NULL else { f <- function(v) { ab <- tryCatch(auc_of(y[idx], score_idx(v$up, v$dn, idx)), error = function(e) NA_real_)
      ao <- tryCatch(auc_of(y, score_idx(v$up, v$dn, ALL)), error = function(e) NA_real_); c(app = ab, orig = ao, thr = v$thr, n = length(v$cf)) }
    list(se = f(m$se), mn = f(m$mn)) } }
bsum <- function(v) { M <- do.call(rbind, lapply(boot, function(z) if (is.null(z)) rep(NA_real_, 4) else z[[v]])); colnames(M) <- c("app", "orig", "thr", "n"); M }
Bse <- bsum("se"); Bmn <- bsum("mn")
# ---- (4) selection frequency: subsamples of 80 % of the samples (stratified by group x batch), whole recipe
sub_idx <- function(s) { set.seed(SUB_SEED0 + s); sort(unlist(lapply(unique(strata_all), function(cl) { ix <- which(strata_all == cl); n <- length(ix)
  ix[sample.int(n, max(1L, min(n, as.integer(round(0.8 * n)))))] }))) }
sel_up <- setNames(numeric(length(feat)), feat); sel_dn <- sel_up; nsub_ok <- 0L
for (s in seq_len(NSUB)) { m <- safe_fit(sub_idx(s)); if (is.null(m)) next; nsub_ok <- nsub_ok + 1L; sel_up[m$se$up] <- sel_up[m$se$up] + 1; sel_dn[m$se$dn] <- sel_dn[m$se$dn] + 1 }
list(ok = TRUE)
}, error = function(e) e)
if (inherits(res_all, "error")) not_run(paste("error in the validation step:", conditionMessage(res_all)))
# ---- summaries
pct <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
row <- function(metric, value, lo, hi, label) data.frame(metric = metric, value = signif(value, 6), ci_low = signif(lo, 6), ci_high = signif(hi, 6), label = label, stringsAsFactors = FALSE)
L_APP <- "apparent (pipeline development scores): OPTIMISTIC, genes and cutoff chosen on these samples"
L_NCV <- sprintf("repeated nested CV (outer %d-fold, %d repeats; LASSO+score+cutoff refitted in every fold; candidate screening NOT repeated): still optimistic; interval = 2.5-97.5 percentile over repeats (repeat-to-repeat variation only, not a confidence interval)", K, NREP)
L_BOOT <- sprintf("Harrell bootstrap (B = %d; LASSO+score+cutoff refitted in every resample; candidate screening NOT repeated): still optimistic", NBOOT)
rows <- list(
  row("pipeline_model_reproduced", as.integer(repro), NA, NA, if (repro) "re-implemented recipe reproduces the pipeline's genes, lambda.1se, cutoff and scores" else "WARNING: re-implemented recipe does NOT reproduce the pipeline model (see log); read the numbers below with caution"),
  row("apparent_AUC", app_auc, NA, NA, L_APP), row("apparent_sensitivity", app_se, NA, NA, L_APP), row("apparent_specificity", app_sp, NA, NA, L_APP))
for (v in c("se", "mn")) {
  d <- if (v == "se") ncv_se else ncv_mn; B <- if (v == "se") Bse else Bmn; tag <- if (v == "se") "" else "_lambda_min"
  lamlab <- if (v == "se") "lambda.1se (the pipeline's rule)" else "SENSITIVITY: lambda.min instead of lambda.1se"
  rows <- c(rows, list(
    row(paste0("nestedCV_AUC", tag), mean(d$auc, na.rm = TRUE), pct(d$auc)[1], pct(d$auc)[2], paste(lamlab, "|", L_NCV)),
    row(paste0("nestedCV_sensitivity", tag), mean(d$sens, na.rm = TRUE), pct(d$sens)[1], pct(d$sens)[2], paste(lamlab, "|", L_NCV)),
    row(paste0("nestedCV_specificity", tag), mean(d$spec, na.rm = TRUE), pct(d$spec)[1], pct(d$spec)[2], paste(lamlab, "|", L_NCV))))
  opt_ <- mean(B[, "app"] - B[, "orig"], na.rm = TRUE)
  apparent_v <- if (v == "se") app_auc else { smn <- score_idx(fm$mn$up, fm$mn$dn, ALL); auc_of(y, smn) }
  rows <- c(rows, list(
    row(paste0("bootstrap_optimism", tag), opt_, NA, NA, paste(lamlab, "| optimism = mean(AUC of bootstrap model on its bootstrap sample - AUC on the original sample) |", L_BOOT)),
    row(paste0("bootstrap_corrected_AUC", tag), apparent_v - opt_, NA, NA, paste(lamlab, "| apparent AUC minus optimism |", L_BOOT)),
    row(paste0("bootstrap_cutoff", tag), mean(B[, "thr"], na.rm = TRUE), pct(B[, "thr"])[1], pct(B[, "thr"])[2], paste(lamlab, "| median-of-folds cutoff of each bootstrap model; value = mean, interval = 2.5-97.5 percentile; the pipeline's own cutoff is", signif(SP_CUT, 8))),
    row(paste0("bootstrap_n_genes", tag), mean(B[, "n"], na.rm = TRUE), pct(B[, "n"])[1], pct(B[, "n"])[2], paste(lamlab, "| number of selected genes per bootstrap model")))) }
rows <- c(rows, list(row("failed_fits", sum(ncv_se$failed) + sum(is.na(Bse[, "app"])) + (NSUB - nsub_ok), NA, NA, "number of resampled fits that errored and were left out")))
S <- do.call(rbind, rows); wtab(S, f_sum)
cand_dir <- setNames(ct$direction_fc, ct$gene_id)
fr <- data.frame(gene_id = feat, symbol = unname(sym[feat]), up = sel_up, dn = sel_dn, stringsAsFactors = FALSE)
fr$frequency <- round((fr$up + fr$dn) / max(nsub_ok, 1), 4); fr$direction <- ifelse(fr$up + fr$dn == 0, NA, ifelse(fr$up >= fr$dn, "up", "down"))
fr$in_final_model <- fr$gene_id %in% gt$gene_id
fr <- fr[fr$frequency > 0 | fr$in_final_model, ]; fr <- fr[order(-fr$frequency, -fr$in_final_model, fr$gene_id), ]
wtab(fr[, c("gene_id", "symbol", "direction", "frequency", "in_final_model")], f_freq)
# ---- INTERNAL_VALIDATION.md
g <- function(m) S[S$metric == m, ]; fmt <- function(m, p = 3) { r <- g(m); if (is.na(r$ci_low)) sprintf("%.*f", p, r$value) else sprintf("%.*f (%.*f to %.*f)", p, r$value, p, r$ci_low, p, r$ci_high) }
fg <- fr[fr$in_final_model, ]
md <- c("# Internal validation of the diagnostic model", "",
 paste0("**", LIMIT, "**"), "",
 sprintf("Samples: %d (%s %d, %s %d). Candidate genes entering the LASSO: %d. The pipeline's own model: %d genes, cutoff %.6g, lambda.1se %.6g.",
         N, base_group, sum(y == 0), pos_group, sum(y == 1), length(feat), nrow(gt), SP_CUT, SP_LAM),
 if (repro) "Check: the re-implemented recipe, run on all samples, reproduces the pipeline's genes, lambda.1se, cutoff and scores." else "**WARNING: the re-implemented recipe did NOT reproduce the pipeline model exactly; read the numbers with caution.**", "",
 "## Numbers (ROC AUC; sensitivity and specificity at the cutoff each model chose)", "",
 "| Estimate | AUC | Sensitivity | Specificity |", "|---|---|---|---|",
 sprintf("| Apparent (pipeline development scores; optimistic) | %s | %s | %s |", fmt("apparent_AUC"), fmt("apparent_sensitivity"), fmt("apparent_specificity")),
 sprintf("| Repeated nested CV, lambda.1se (%d repeats; interval = 2.5-97.5 percentile over repeats) | %s | %s | %s |", NREP, fmt("nestedCV_AUC"), fmt("nestedCV_sensitivity"), fmt("nestedCV_specificity")),
 sprintf("| Repeated nested CV, lambda.min (sensitivity analysis) | %s | %s | %s |", fmt("nestedCV_AUC_lambda_min"), fmt("nestedCV_sensitivity_lambda_min"), fmt("nestedCV_specificity_lambda_min")), "",
 sprintf("Bootstrap (B = %d, Harrell): optimism %s, optimism-corrected AUC %s (lambda.min: optimism %s, corrected AUC %s).",
         NBOOT, fmt("bootstrap_optimism"), fmt("bootstrap_corrected_AUC"), fmt("bootstrap_optimism_lambda_min"), fmt("bootstrap_corrected_AUC_lambda_min")),
 sprintf("Cutoff over bootstrap models: mean %s (2.5-97.5 percentile interval shown), the pipeline's cutoff is %.6g. Genes selected per bootstrap model: %s.", fmt("bootstrap_cutoff"), SP_CUT, fmt("bootstrap_n_genes", 1)), "",
 sprintf("## Selection frequency (%d subsamples of 80 %% of the samples, stratified; fraction of subsamples in which the gene is selected at lambda.1se)", NSUB), "",
 "Report only: no consensus panel is built and the model is not changed.", "",
 "| Gene | Direction | Frequency | In final model |", "|---|---|---|---|",
 sprintf("| %s | %s | %.2f | %s |", ifelse(is.na(fr$symbol) | fr$symbol == "", fr$gene_id, fr$symbol), fr$direction, fr$frequency, ifelse(fr$in_final_model, "yes", "no"))[seq_len(min(nrow(fr), 25))], "",
 if (nrow(fr) > 25) sprintf("(first 25 of %d genes; the full table is selection_frequency.tsv)", nrow(fr)) else NULL, "",
 "## Method", "",
 "The whole model-building recipe of the pipeline (fit_diagnostic_model.R) is repeated inside every resample: log2(CPM + 1); glmnet binomial, alpha 1, folds stratified by group x batch, lambda.1se; up and down gene sets from the coefficient signs; score = mean(up) - mean(down); cutoff = median of the inner-fold Youden thresholds. The candidate list is fixed (the genes that passed rules 1-3b on all samples).",
 sprintf("Nested CV: outer %d-fold stratified by group x batch, %d repeats, out-of-fold scores pooled per repeat for the AUC; each outer fold uses its own model and cutoff. Bootstrap: %d resamples drawn within each class; optimism = mean(AUC of the bootstrap model on its bootstrap sample - AUC on the original sample); corrected AUC = apparent AUC - optimism. The bootstrap models are scored on the original samples that include their own training samples, as in Harrell's procedure.", K, NREP, NBOOT),
 sprintf("Seeds (fixed): the model's own seeds %d (LASSO folds) and %d (inner cutoff folds) are reused inside every fit; outer CV repeat r uses %d + r, bootstrap b uses %d + b, subsample s uses %d + s.", SEED, CSEED, OUTER_SEED0, BOOT_SEED0, SUB_SEED0),
 sprintf("Versions: R %s, glmnet %s, pROC %s.", getRversion(), packageVersion("glmnet"), packageVersion("pROC")), "",
 "## What this does not show", "",
 LIMIT, "Small samples give wide, unstable estimates; with the same batch and patient structure in every resample, batch effects are not tested here. Only an external cohort measured separately can validate the biomarker.")
writeLines(md[!is.null(md)], f_md)
cat("internal validation written to", normalizePath(OUT), "\n")
