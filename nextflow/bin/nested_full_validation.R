#!/usr/bin/env Rscript
# Fully NESTED internal validation of the diagnostic signature (candidate screen INSIDE every training fold). ADDITIVE and REPORT-ONLY: it never feeds
# back into the model and changes no existing file or output. Companion of internal_validation.R, which repeats only LASSO + score + cutoff and keeps the
# candidate screen outside the resamples (optimistic). Here, in every outer training fold and using the TRAINING samples only:
#   1. DESeq2 (the deseq2.R recipe: ~ batch + group, or ~ group with one batch; pre-filter and apeglm shrinkage from conf/thresholds.config)
#   2. candidate rules 1 and 2 exactly as candidate_rules.R applies them (padj < padj_max, |apeglm-shrunk log2FC| >= abs_shrunk_log2fc_min, baseMean >= basemean_min)
#   3. LASSO (glmnet binomial, alpha 1, lambda.1se), up/down score, cutoff = median of K inner-fold Youden cutoffs (fit_diagnostic_model.R recipe)
#   4. the HELD-OUT fold is scored with the training-fold model and the training-fold cutoff.
# NOT RUN here: rule 3a (needs the second quantification route) and rule 3b (leave-one-out: N more DESeq2 fits per fold). The nested estimate therefore
# describes the SIMPLIFIED procedure "rules 1-2 + LASSO + score + cutoff", not the full locked procedure. It is NOT independent validation.
# A label-permutation null (labels permuted within batch, the ENTIRE procedure repeated) gives the null AUC distribution and an empirical p-value.
# Packages: DESeq2 + apeglm live in the r_dge environment, glmnet + pROC in r_classifier; no environment has all four. The script therefore
# orchestrates itself in two stages (--stage screen under r_dge, --stage model under r_classifier); the default --stage all needs only base R.
# usage: nested_full_validation.R --counts counts_training.tsv --sheet sample_sheet.discover.tsv --questions config/questions.yaml --rules config/candidate_rules.yaml
#        [--thresholds conf/thresholds.config] [--outer-folds 5] [--repeats 5] [--permutations 20] [--perm-repeats <= repeats> ] [--cores 2 (max 4)] [--seed 1000]
#        [--out .] [--patient-level-cv true] [--check-only true] [--keep-work false] [--rscript-dge <path>] [--rscript-clf <path>]
#        [--screen-outside true]  COMPARATOR ONLY: screen once on ALL samples, then nested LASSO+cutoff (the internal_validation.R design); permutations are not run.
# Exit status: 0 always for analysis problems (a 'NOT RUN: reason' row is written); 2 for usage errors; --check-only exits 1 when the inputs/environment are not usable.
a <- commandArgs(trailingOnly = TRUE)
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
KNOWN <- c("counts", "sheet", "questions", "rules", "thresholds", "outer-folds", "repeats", "permutations", "perm-repeats", "cores", "seed", "out", "patient-level-cv",
           "check-only", "keep-work", "rscript-dge", "rscript-clf", "stage", "screen-outside")
usage_error <- function(msg) { cat("USAGE ERROR:", msg, "\n", file = stderr()); quit(save = "no", status = 2) }
opt <- list()
if (length(a) %% 2 != 0) usage_error("options must be given as '--name value' pairs")
for (i in seq(1, length(a), by = 2)) { if (!grepl("^--", a[i])) usage_error(paste("unexpected argument", a[i])); nm <- sub("^--", "", a[i]); if (!nm %in% KNOWN) usage_error(paste("unknown option --", nm, sep = "")); opt[[nm]] <- a[i + 1] }
for (k in c("counts", "sheet", "questions", "rules")) if (is.null(opt[[k]])) usage_error(paste0("missing --", k))
int_opt <- function(nm, def, lo) { v <- if (is.null(opt[[nm]])) def else suppressWarnings(as.integer(opt[[nm]])); if (is.na(v) || v < lo) usage_error(sprintf("--%s must be an integer >= %d", nm, lo)); v }
KO <- int_opt("outer-folds", 5L, 2L); NREP <- int_opt("repeats", 5L, 1L); NPERM <- int_opt("permutations", 20L, 0L); SEED <- int_opt("seed", 1000L, 0L)
PREP <- int_opt("perm-repeats", NREP, 1L); CORES_REQ <- int_opt("cores", 2L, 1L); MAXC <- suppressWarnings(as.integer(Sys.getenv("NESTED_MAX_CORES", "4"))); if (is.na(MAXC) || MAXC < 1L) MAXC <- 4L; CORES <- min(CORES_REQ, MAXC)   # cap 4 on a shared pod; raise it with NESTED_MAX_CORES on an idle one
bool_opt <- function(nm) { v <- tolower(if (is.null(opt[[nm]])) "false" else opt[[nm]]); if (!v %in% c("true", "false")) usage_error(paste0("--", nm, " must be true or false")); v == "true" }
PLC <- if (is.null(opt[["patient-level-cv"]])) TRUE else bool_opt("patient-level-cv"); KEEP <- bool_opt("keep-work"); OUTSIDE <- bool_opt("screen-outside"); CHECK <- bool_opt("check-only")
STAGE <- if (is.null(opt$stage)) "all" else opt$stage; if (!STAGE %in% c("all", "screen", "model")) usage_error("--stage must be all, screen or model")
if (OUTSIDE) NPERM <- 0L
OUT <- if (is.null(opt$out)) "." else opt$out
RS_DGE <- if (is.null(opt[["rscript-dge"]])) "/workspace/Biomarker_pipeline/env/prefix/r_dge/bin/Rscript" else opt[["rscript-dge"]]
RS_CLF <- if (is.null(opt[["rscript-clf"]])) "/workspace/Biomarker_pipeline/env/prefix/r_classifier/bin/Rscript" else opt[["rscript-clf"]]
THR_FILE <- if (!is.null(opt$thresholds)) opt$thresholds else file.path(dirname(dirname(dirname(normalizePath(SELF, mustWork = FALSE)))), "conf", "thresholds.config")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
.kf <- tempfile(); writeLines(c(paste(unname(tools::md5sum(c(opt$counts, opt$sheet, opt$questions, opt$rules, THR_FILE))), collapse = ","), KO, NREP, NPERM, PREP, SEED, PLC, OUTSIDE), .kf)
WORK <- file.path(OUT, paste0("nested_work_", substr(unname(tools::md5sum(.kf)), 1, 8)))      # resume cache, keyed on the inputs and settings; removed after a successful run unless --keep-work true
TAG <- "NESTED, rules 1-2 only; NOT independent validation"
if (OUTSIDE) TAG <- "COMPARATOR, screen on ALL samples (rules 1-2) OUTSIDE the resampling; optimistic BY DESIGN; NOT nested, NOT independent validation"
LIMIT <- paste("LIMITATION: this nested estimate repeats candidate rules 1 and 2, the LASSO, the score and the cutoff inside every training fold, but rule 3a (second quantification route) and",
               "rule 3b (leave-one-out robustness) are NOT run, so it estimates the SIMPLIFIED procedure 'rules 1-2 + LASSO + score + cutoff', not the full locked procedure.",
               "Pooled out-of-fold scores come from different fold models with different cutoffs. With few samples every number is wide and unstable.",
               "This is internal validation, NOT independent validation; only an external cohort measured separately can validate the signature.")
f_sum <- file.path(OUT, "nested_validation_summary.tsv"); f_freq <- file.path(OUT, "nested_selection_frequency.tsv")
f_null <- file.path(OUT, "nested_permutation_null.tsv"); f_md <- file.path(OUT, "NESTED_VALIDATION.md")
wtab <- function(d, f) write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
not_run <- function(reason) {
  lab <- paste0("NOT RUN: ", reason)
  wtab(data.frame(metric = "nested_validation", value = NA, ci_low = NA, ci_high = NA, label = lab), f_sum)
  wtab(data.frame(gene_id = lab, n_selected = NA, n_folds_run = NA, frequency = NA, candidate_frequency = NA, direction = NA), f_freq)
  wtab(data.frame(permutation = NA, mean_pooled_auc = NA, folds_run = NA, folds_skipped = NA, note = lab), f_null)
  writeLines(c("# Fully nested internal validation", "", paste0("**", lab, "**"), "",
               "No nested estimate was produced. The diagnostic model and every other pipeline output are unaffected.", "", LIMIT), f_md)
  cat(lab, "\n"); quit(save = "no", status = 0)
}
pkg_ok <- function(p) requireNamespace(p, quietly = TRUE)
# ---------------------------------------------------------------- settings read like the pipeline scripts do
yaml_vals <- function(file, block, key) { l <- readLines(file); i <- grep(paste0("^", block, ":"), l); stopifnot(length(i) == 1)
  j <- i + 1; while (j <= length(l) && grepl("^(\\s|#|$)", l[j])) j <- j + 1; blk <- l[(i + 1):(j - 1)]
  v <- grep(paste0("^\\s*", key, ":"), blk, value = TRUE); as.numeric(sub("^\\s*[a-z0-9_]+:\\s*([0-9.]+).*$", "\\1", v)) }
read_key <- function(file, key, sep) { l <- grep(paste0("^\\s*", key, "\\s*", sep), readLines(file), value = TRUE); stopifnot(length(l) == 1)
  trimws(gsub("['\"]", "", sub("(#|//).*$", "", sub(paste0("^\\s*", key, "\\s*", sep), "", l)))) }
# ---------------------------------------------------------------- inputs and the plan (base R only: identical in both stages)
load_inputs <- function() {
  for (f in c(opt$counts, opt$sheet, opt$questions, opt$rules, THR_FILE)) if (!file.exists(f)) stop("file not found: ", f)
  S <- list()
  S$P_CUT <- yaml_vals(opt$rules, "rule1_significance", "padj_max"); S$LFC_CUT <- yaml_vals(opt$rules, "rule1_significance", "abs_shrunk_log2fc_min")
  bm <- yaml_vals(opt$rules, "rule2_measurable", "basemean_min"); S$BM_MIN <- bm[1]
  stopifnot(length(S$P_CUT) == 1, length(S$LFC_CUT) == 1, is.finite(c(S$P_CUT, S$LFC_CUT, S$BM_MIN)))
  S$ALPHA <- as.numeric(read_key(THR_FILE, "de_padj", "=")); S$MIN_COUNT <- as.numeric(read_key(THR_FILE, "de_min_count", "=")); S$DE_MIN_SAMPLES <- read_key(THR_FILE, "de_min_samples", "=")
  S$LSEED <- as.integer(read_key(opt$questions, "lasso_seed", ":")); S$CSEED <- as.integer(read_key(opt$questions, "cutoff_seed", ":")); S$KIN <- as.integer(read_key(opt$questions, "cv_folds", ":"))
  stopifnot(is.finite(c(S$ALPHA, S$MIN_COUNT, S$LSEED, S$CSEED, S$KIN)))
  cn <- as.matrix(read.delim(opt$counts, row.names = 1, check.names = FALSE))
  ss <- read.delim(opt$sheet, stringsAsFactors = FALSE, colClasses = "character"); rownames(ss) <- ss$sample_id
  if (!all(colnames(cn) %in% ss$sample_id)) stop("samples of the matrix missing from the sample sheet")
  ss <- ss[colnames(cn), ]
  base_group <- unique(ss$group[toupper(ss$baseline) == "TRUE"]); if (length(base_group) != 1) stop("need exactly one baseline group")
  pos_group <- setdiff(unique(ss$group), base_group); if (length(pos_group) != 1) stop("need exactly two groups")
  S$cn <- cn; S$base <- base_group; S$pos <- pos_group; S$y <- setNames(as.integer(ss$group == pos_group), colnames(cn)); S$batch <- setNames(ss$batch, colnames(cn))
  S$pat <- setNames(if ("patient_id" %in% names(ss)) ss$patient_id else colnames(cn), colnames(cn)); S$N <- ncol(cn)
  S$REP_PAT <- PLC && any(duplicated(S$pat)); S
}
strat_foldid <- function(strata, k) { id <- integer(length(strata)); for (cl in unique(strata)) { ix <- which(strata == cl); id[ix] <- sample(rep_len(seq_len(k), length(ix))) }; id }
strat_foldid_pat <- function(strata, patient, k) {
  up <- unique(patient); sig <- vapply(up, function(p) paste(sort(unique(strata[patient == p])), collapse = "|"), "")
  pid <- integer(length(up)); names(pid) <- up
  for (cl in unique(sig)) { ix <- which(sig == cl); pid[ix] <- sample(rep_len(seq_len(k), length(ix))) }
  unname(pid[as.character(patient)]) }
make_foldid <- function(S, strata, patient, k) if (S$REP_PAT) strat_foldid_pat(strata, patient, k) else strat_foldid(strata, k)
permute_labels <- function(S, b) {      # labels permuted WITHIN batch (patient-level when a patient_id repeats and each patient has one group and one batch)
  set.seed(S$seed0 + 7919L * b); y <- S$y
  pl <- tapply(paste(y, S$batch), S$pat, function(v) length(unique(v)) == 1)
  if (S$REP_PAT && all(pl)) {
    ps <- !duplicated(S$pat); pb <- S$batch[ps]; py <- setNames(y[ps], S$pat[ps]); pb <- setNames(pb, S$pat[ps])
    for (bt in unique(pb)) { ix <- which(pb == bt); py[ix] <- py[ix][sample.int(length(ix))] }
    return(setNames(unname(py[as.character(S$pat)]), names(y))) }
  for (bt in unique(S$batch)) { ix <- which(S$batch == bt); y[ix] <- y[ix][sample.int(length(ix))] }
  y }
build_plan <- function(S) {
  S$seed0 <- SEED; labs <- list(); tasks <- list()
  for (b in 0:NPERM) {
    yb <- if (b == 0) S$y else permute_labels(S, b); labs[[as.character(b)]] <- yb; nr <- if (b == 0) NREP else PREP
    for (r in seq_len(nr)) {
      set.seed(SEED + 1000L * b + r); ofid <- make_foldid(S, paste(yb, S$batch), S$pat, KO)
      for (k in seq_len(KO)) {
        te <- which(ofid == k); tr <- which(ofid != k); why <- ""
        tb <- table(S$batch[tr], yb[tr]); tcl <- table(factor(yb[tr], 0:1))
        if (length(te) == 0 || length(tr) == 0) why <- "empty fold"
        else if (any(tcl < 2)) why <- "a class has fewer than 2 training samples"
        else if (ncol(tb) < 2 || any(tb == 0)) why <- "training fold has an empty batch x group cell"
        tasks[[length(tasks) + 1]] <- list(id = sprintf("p%d_r%d_f%d", b, r, k), b = b, r = r, k = k, tr = tr, te = te, skip = why) } } }
  list(labs = labs, tasks = tasks) }
scr_file <- function(id) file.path(WORK, paste0("screen_", id, ".rds"))
# ---------------------------------------------------------------- STAGE screen (r_dge): DESeq2 + rules 1-2 on training samples
stage_screen <- function() {
  suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(parallel) })
  S <- load_inputs(); pl <- build_plan(S); cn <- S$cn; dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
  one <- function(cols, yv, id) {
    f <- scr_file(id); if (file.exists(f)) return(invisible(NULL))
    res <- tryCatch({
      cd <- data.frame(group = factor(ifelse(yv[cols] == 1, S$pos, S$base), levels = c(S$base, S$pos)), batch = factor(S$batch[cols]), row.names = colnames(cn)[cols])
      design <- if (nlevels(cd$batch) > 1) ~ batch + group else ~ group
      stopifnot(all(table(cd$batch, cd$group) > 0))
      d <- DESeqDataSetFromMatrix(round(cn[, cols, drop = FALSE]), colData = cd, design = design)
      ms <- if (S$DE_MIN_SAMPLES == "auto") min(table(cd$group)) else as.numeric(S$DE_MIN_SAMPLES)      # 'auto' = smallest group of the TRAINING samples, as deseq2.R does on its input
      d <- d[rowSums(counts(d) >= S$MIN_COUNT) >= ms, ]; d <- DESeq(d, quiet = TRUE)
      coef <- paste0("group_", S$pos, "_vs_", S$base); stopifnot(coef %in% resultsNames(d))
      rs <- results(d, name = coef, alpha = S$ALPHA); sh <- lfcShrink(d, coef = coef, type = "apeglm", res = rs, quiet = TRUE)
      r1 <- !is.na(rs$padj) & rs$padj < S$P_CUT & !is.na(sh$log2FoldChange) & abs(sh$log2FoldChange) >= S$LFC_CUT      # rule 1 (candidate_rules.R r1)
      r2 <- !is.na(rs$baseMean) & rs$baseMean >= S$BM_MIN                                                             # rule 2, single route
      list(ok = TRUE, n_tested = sum(!is.na(rs$padj)), n_rule1 = sum(r1), cand = data.frame(gene_id = rownames(rs)[r1 & r2], dir = sign(sh$log2FoldChange[r1 & r2]), stringsAsFactors = FALSE))
    }, error = function(e) list(ok = FALSE, reason = paste("DESeq2 screen failed:", conditionMessage(e))))
    saveRDS(res, f) }
  if (OUTSIDE) { cat("screen on ALL samples (comparator)\n"); one(seq_len(S$N), pl$labs[["0"]], "all"); return(invisible(NULL)) }
  todo <- Filter(function(t) t$skip == "", pl$tasks); cat(sprintf("screen: %d DESeq2 fits (%d tasks skipped by the plan), %d cores\n", length(todo), length(pl$tasks) - length(todo), CORES))
  invisible(mclapply(todo, function(t) { one(t$tr, pl$labs[[as.character(t$b)]], t$id); NULL }, mc.cores = CORES, mc.preschedule = FALSE)) }
# ---------------------------------------------------------------- STAGE model (r_classifier): LASSO + score + cutoff + held-out prediction + summaries
stage_model <- function() {
  suppressPackageStartupMessages({ library(glmnet); library(pROC); library(parallel) }); options(warn = -1)
  S <- load_inputs(); pl <- build_plan(S); cn <- S$cn; y0 <- S$y
  L4 <- log2(sweep(cn, 2, colSums(cn), "/") * 1e6 + 1)      # log2(CPM + 1), CPM total = column sum (a per-sample quantity: no leakage between samples)
  if (S$N < 10) not_run(sprintf("only %d samples (at least 10 needed)", S$N))
  if (min(table(y0)) < 5) not_run(sprintf("a class has fewer than 5 samples (%s %d, %s %d)", S$base, sum(y0 == 0), S$pos, sum(y0 == 1)))
  youden_thr <- function(yy, s) { r <- roc(yy, s, levels = c(0, 1), direction = "<", quiet = TRUE); th <- coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE)[, 1]; median(th[is.finite(th)]) }
  score_fun <- function(up, dn, idx) vapply(idx, function(j) { u <- if (length(up)) mean(L4[up, j]) else 0; d <- if (length(dn)) mean(L4[dn, j]) else 0; u - d }, 0)
  cand_all <- if (OUTSIDE) readRDS(scr_file("all")) else NULL
  run_task <- function(t) {
    yb <- pl$labs[[as.character(t$b)]]; base <- list(id = t$id, b = t$b, r = t$r, k = t$k, te = t$te, status = "skipped", reason = t$skip, n_cand = NA_integer_, n_sel = NA_integer_)
    if (t$skip != "") return(base)
    sc <- if (OUTSIDE) cand_all else { f <- scr_file(t$id); if (file.exists(f)) readRDS(f) else list(ok = FALSE, reason = "screen result missing") }
    if (!isTRUE(sc$ok)) { base$reason <- sc$reason; return(base) }
    feat <- intersect(sc$cand$gene_id, rownames(L4)); base$n_cand <- length(feat)
    if (length(feat) < 2) { base$reason <- sprintf("fewer than 2 candidates after rules 1-2 (%d)", length(feat)); return(base) }
    tryCatch({
      tr <- t$tr; yt <- yb[tr]; strata_tr <- paste(yt, S$batch[tr]); patt <- S$pat[tr]
      X <- t(L4[feat, tr, drop = FALSE])
      set.seed(S$LSEED); fid <- make_foldid(S, strata_tr, patt, S$KIN)
      cv <- cv.glmnet(X, yt, family = "binomial", alpha = 1, foldid = fid, type.measure = "deviance")
      cf <- as.matrix(coef(cv, s = "lambda.1se")); cf <- cf[rownames(cf) != "(Intercept)", 1]; cf <- cf[cf != 0]
      base$n_sel <- length(cf)
      if (length(cf) == 0) { base$reason <- "LASSO selected no gene at lambda.1se"; return(base) }
      up <- names(cf)[cf > 0]; dn <- names(cf)[cf < 0]
      set.seed(S$CSEED); fid2 <- make_foldid(S, strata_tr, patt, S$KIN); s_tr <- score_fun(up, dn, tr); thr_f <- numeric(S$KIN)
      for (kk in seq_len(S$KIN)) { itr <- fid2 != kk; thr_f[kk] <- youden_thr(yt[itr], s_tr[itr]) }
      thr <- median(thr_f); s_te <- score_fun(up, dn, t$te)
      base$status <- "ok"; base$reason <- ""; base$score <- s_te; base$thr <- thr; base$genes <- cf; base$cand_genes <- sc$cand; base }, 
      error = function(e) { base$reason <- paste("model fit failed:", conditionMessage(e)); base }) }
  cat(sprintf("model: %d tasks, %d cores\n", length(pl$tasks), CORES))
  R <- mclapply(pl$tasks, run_task, mc.cores = CORES, mc.preschedule = TRUE)
  bad <- vapply(R, function(z) !is.list(z), TRUE); if (any(bad)) stop("a model worker failed: ", paste(unlist(R[bad]), collapse = "; "))
  # ---- per (permutation b, repeat r): pooled out-of-fold AUC (+ DeLong CI), sensitivity / specificity at each fold's own cutoff (+ Clopper-Pearson)
  cp <- function(x, n) if (n > 0) { ci <- binom.test(x, n)$conf.int; c(x / n, ci[1], ci[2]) } else c(NA, NA, NA)
  pr_rows <- list()
  for (b in 0:NPERM) { yb <- pl$labs[[as.character(b)]]; nr <- if (b == 0) NREP else PREP
    for (r in seq_len(nr)) {
      Z <- Filter(function(z) z$b == b && z$r == r, R); ok <- Filter(function(z) z$status == "ok", Z)
      sc <- rep(NA_real_, S$N); pd <- rep(NA_integer_, S$N)
      for (z in ok) { sc[z$te] <- z$score; pd[z$te] <- as.integer(z$score > z$thr) }
      has <- !is.na(sc); row <- data.frame(b = b, r = r, folds_run = length(ok), folds_skipped = length(Z) - length(ok), n_scored = sum(has),
        auc = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_, TP = NA, FN = NA, TN = NA, FP = NA, mean_genes = if (length(ok)) mean(vapply(ok, function(z) z$n_sel, 0)) else NA_real_,
        mean_cand = if (length(ok)) mean(vapply(ok, function(z) z$n_cand, 0)) else NA_real_)
      if (sum(has) > 1 && length(unique(yb[has])) == 2) {
        ro <- roc(yb[has], sc[has], levels = c(0, 1), direction = "<", quiet = TRUE); ci <- suppressWarnings(ci.auc(ro, method = "delong"))
        row$auc <- as.numeric(auc(ro)); row$auc_lo <- ci[1]; row$auc_hi <- ci[3]
        row$TP <- sum(pd[has] == 1 & yb[has] == 1); row$FN <- sum(pd[has] == 0 & yb[has] == 1); row$TN <- sum(pd[has] == 0 & yb[has] == 0); row$FP <- sum(pd[has] == 1 & yb[has] == 0) }
      pr_rows[[length(pr_rows) + 1]] <- row } }
  PR <- do.call(rbind, pr_rows); O <- PR[PR$b == 0, ]; NL <- PR[PR$b > 0, ]
  # ---- permutation null (one value per permutation = mean over its repeats of the pooled AUC; same statistic as the observed value)
  nullt <- do.call(rbind, lapply(seq_len(NPERM), function(b) { d <- NL[NL$b == b, ]
    data.frame(permutation = b, mean_pooled_auc = if (all(is.na(d$auc))) NA_real_ else mean(d$auc, na.rm = TRUE), n_repeats_with_auc = sum(!is.na(d$auc)), folds_run = sum(d$folds_run), folds_skipped = sum(d$folds_skipped),
               mean_genes = if (all(is.na(d$mean_genes))) NA_real_ else mean(d$mean_genes, na.rm = TRUE)) }))
  if (is.null(nullt)) nullt <- data.frame(permutation = integer(), mean_pooled_auc = numeric(), n_repeats_with_auc = integer(), folds_run = integer(), folds_skipped = integer(), mean_genes = numeric())
  obs_auc <- if (all(is.na(O$auc))) NA_real_ else mean(O$auc, na.rm = TRUE)
  nv <- nullt$mean_pooled_auc[!is.na(nullt$mean_pooled_auc)]; B <- length(nv)
  pval <- if (B > 0 && !is.na(obs_auc)) (sum(nv >= obs_auc) + 1) / (B + 1) else NA_real_
  # ---- selection frequency over the observed (b = 0) folds that ran
  Zobs <- Filter(function(z) z$b == 0 && z$status == "ok", R); nrun <- length(Zobs)
  genes <- unique(unlist(lapply(Zobs, function(z) names(z$genes))))
  fr <- data.frame(gene_id = genes, stringsAsFactors = FALSE)
  fr$n_selected <- as.numeric(vapply(genes, function(g) sum(vapply(Zobs, function(z) g %in% names(z$genes), TRUE)), 0))
  fr$n_folds_run <- rep(nrun, nrow(fr)); fr$frequency <- round(fr$n_selected / max(nrun, 1), 4)
  fr$candidate_frequency <- round(as.numeric(vapply(genes, function(g) sum(vapply(Zobs, function(z) g %in% z$cand_genes$gene_id, TRUE)), 0)) / max(nrun, 1), 4)
  fr$direction <- as.character(vapply(genes, function(g) { s <- unlist(lapply(Zobs, function(z) if (g %in% names(z$genes)) sign(z$genes[[g]]))); if (sum(s > 0) >= sum(s < 0)) "up" else "down" }, ""))
  fr <- fr[order(-fr$frequency, fr$gene_id), ]; rownames(fr) <- NULL
  # ---- summary table
  fl <- function(x) signif(x, 6); row <- function(m, v, lo, hi, lab) data.frame(metric = m, value = fl(v), ci_low = fl(lo), ci_high = fl(hi), label = paste0(lab, " | ", TAG), stringsAsFactors = FALSE)
  nm <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  se_ci <- t(mapply(function(x, n) cp(x, n), O$TP, O$TP + O$FN)); sp_ci <- t(mapply(function(x, n) cp(x, n), O$TN, O$TN + O$FP))
  Lauc <- sprintf("pooled out-of-fold AUC per repeat, then mean over %d repeats; interval = mean of the per-repeat DeLong 95%% CI bounds (DeLong treats the pooled scores as one model: approximate)", NREP)
  Lsp <- "sensitivity/specificity of the pooled out-of-fold calls, each fold at its own training-derived cutoff; mean over repeats; interval = mean of the per-repeat Clopper-Pearson 95% bounds (repeats reuse the same samples: not independent)"
  S_ <- list(
    row("procedure", NA, NA, NA, "rules 1-2 (DESeq2 padj, |apeglm LFC|, baseMean) + LASSO lambda.1se + up/down score + median inner-fold Youden cutoff; rule 3a and rule 3b NOT run: estimates the SIMPLIFIED procedure, not the full locked procedure"),
    row("nested_AUC", obs_auc, nm(O$auc_lo), nm(O$auc_hi), Lauc),
    row("nested_AUC_repeat_spread", obs_auc, suppressWarnings(min(O$auc, na.rm = TRUE)), suppressWarnings(max(O$auc, na.rm = TRUE)), "value = mean over repeats; interval = min and max of the per-repeat pooled AUCs (repeat-to-repeat variation only, NOT a confidence interval)"),
    row("nested_sensitivity", nm(O$TP / (O$TP + O$FN)), nm(se_ci[, 2]), nm(se_ci[, 3]), Lsp),
    row("nested_specificity", nm(O$TN / (O$TN + O$FP)), nm(sp_ci[, 2]), nm(sp_ci[, 3]), Lsp),
    row("folds_run", sum(O$folds_run), NA, NA, sprintf("outer folds that ran, summed over %d repeats x %d folds", NREP, KO)),
    row("folds_skipped", sum(O$folds_skipped), NA, NA, "outer folds skipped (empty batch x group cell in the training fold, fewer than 2 candidates, no LASSO gene, or a failed fit); reasons in NESTED_VALIDATION.md"),
    row("mean_selected_genes", nm(O$mean_genes), NA, NA, "mean number of genes selected by the LASSO per fold that ran"),
    row("mean_candidates_after_rules_1_2", nm(O$mean_cand), NA, NA, "mean number of candidates entering the LASSO per fold that ran"))
  for (i in seq_len(nrow(O))) S_[[length(S_) + 1]] <- row(sprintf("nested_AUC_repeat_%d", O$r[i]), O$auc[i], O$auc_lo[i], O$auc_hi[i], "pooled out-of-fold AUC of this repeat with DeLong 95% CI")
  if (!OUTSIDE) {
    S_ <- c(S_, list(row("null_AUC_mean", nm(nv), NA, NA, sprintf("label-permutation null (B = %d usable of %d permutations; labels permuted within batch; whole procedure repeated, %d repeats each)", B, NPERM, PREP)),
      row("null_AUC_p95", if (B) quantile(nv, 0.95, names = FALSE) else NA, NA, NA, "95th percentile of the null AUC distribution"),
      row("permutation_p_value", pval, NA, NA, sprintf("empirical p = (b+1)/(B+1), b = null AUCs >= observed %.4f, B = %d; the smallest possible value is 1/(B+1)", obs_auc, B))))
  }
  S_ <- do.call(rbind, S_); if (!OUTSIDE) wtab(nullt, f_null) else wtab(data.frame(note = "comparator mode (--screen-outside true): permutations not run"), f_null)
  wtab(S_, f_sum); wtab(fr, f_freq)
  # ---- NESTED_VALIDATION.md
  Zsk <- Filter(function(z) z$b == 0 && z$status != "ok", R); rs <- if (length(Zsk)) sort(table(vapply(Zsk, function(z) z$reason, "")), decreasing = TRUE) else NULL
  g <- function(m) S_[S_$metric == m, ]; f3 <- function(m) { x <- g(m); if (is.na(x$value)) "NA" else if (is.na(x$ci_low)) sprintf("%.3f", x$value) else sprintf("%.3f (%.3f to %.3f)", x$value, x$ci_low, x$ci_high) }
  md <- c(if (OUTSIDE) "# COMPARATOR: screen on ALL samples, then nested LASSO (optimistic by design)" else "# Fully nested internal validation of the diagnostic signature", "",
    paste0("**", LIMIT, "**"), "",
    if (OUTSIDE) c("**This is the COMPARATOR mode: the candidate screen (rules 1-2) was run ONCE on all samples and only LASSO + score + cutoff were repeated inside the folds, as internal_validation.R does. The numbers are optimistic by design and must not be reported as a validation result.**", "") else NULL,
    sprintf("Samples: %d (%s %d, %s %d). Outer %d-fold cross-validation, %d repeats; %d label permutations (%d repeats each). Fold assignment stratified by group x batch%s. Seed %d.",
            S$N, S$base, sum(y0 == 0), S$pos, sum(y0 == 1), KO, NREP, NPERM, PREP, if (S$REP_PAT) ", all samples of one patient kept together" else "", SEED), "",
    "## What was done in every training fold", "",
    "1. DESeq2 on the training samples only (design ~ batch + group, or ~ group with one batch; pre-filter and apeglm shrinkage as in the main pipeline).",
    sprintf("2. Candidate rules 1 and 2: padj < %g, |apeglm-shrunk log2FC| >= %g, baseMean >= %g. Rule 3a and rule 3b are NOT run (cost and the second quantification route).", S$P_CUT, S$LFC_CUT, S$BM_MIN),
    "3. LASSO (glmnet, binomial, lambda.1se), up/down score (mean of up genes minus mean of down genes, log2 CPM + 1), cutoff = median of the inner-fold Youden thresholds.",
    "4. The held-out fold is scored with the training-fold model and the training-fold cutoff. Out-of-fold scores are pooled per repeat.", "",
    "## Numbers", "", "| Estimate | Value |", "|---|---|",
    sprintf("| Nested AUC, mean over repeats (DeLong CI of the pooled scores) | %s |", f3("nested_AUC")),
    sprintf("| Spread of the per-repeat AUC (min to max) | %s |", f3("nested_AUC_repeat_spread")),
    sprintf("| Sensitivity (Clopper-Pearson) | %s |", f3("nested_sensitivity")), sprintf("| Specificity (Clopper-Pearson) | %s |", f3("nested_specificity")),
    sprintf("| Folds run / skipped (observed, all repeats) | %d / %d |", sum(O$folds_run), sum(O$folds_skipped)),
    sprintf("| Mean genes selected by the LASSO per fold | %s |", if (is.na(nm(O$mean_genes))) "NA" else sprintf("%.1f", nm(O$mean_genes))),
    sprintf("| Mean candidates after rules 1-2 per fold | %s |", if (is.na(nm(O$mean_cand))) "NA" else sprintf("%.1f", nm(O$mean_cand))), "",
    if (!OUTSIDE) c("## Permutation null", "",
      sprintf("Labels were permuted within batch and the whole procedure was repeated %d times (usable: %d). Null AUC: mean %s, 95th percentile %s. Observed AUC %s gives an empirical p-value of %s = (b+1)/(B+1); with B = %d the smallest possible p-value is %.3f.",
              NPERM, B, if (B) sprintf("%.3f", mean(nv)) else "NA", if (B) sprintf("%.3f", quantile(nv, 0.95, names = FALSE)) else "NA", if (is.na(obs_auc)) "NA" else sprintf("%.3f", obs_auc), if (is.na(pval)) "NA" else sprintf("%.3f", pval), B, 1 / (B + 1)),
      "A null mean clearly above 0.5 would signal residual selection bias or fold-structure leakage in the procedure.", "") else NULL,
    "## Skipped folds (observed run)", "", if (is.null(rs)) "No fold was skipped." else sprintf("- %d x %s", as.integer(rs), names(rs)), "",
    if (!OUTSIDE && NPERM > 0) { Zp <- Filter(function(z) z$b > 0 && z$status != "ok", R); rp <- if (length(Zp)) sort(table(vapply(Zp, function(z) z$reason, "")), decreasing = TRUE) else NULL
      c("## Skipped folds (permutation runs)", "", if (is.null(rp)) "No fold was skipped." else sprintf("- %d x %s", as.integer(rp), names(rp)), "") } else NULL,
    sprintf("## Selection frequency (fraction of the %d observed folds that ran in which the LASSO selected the gene; top %d, full table nested_selection_frequency.tsv)", nrun, min(nrow(fr), 20)), "",
    if (nrow(fr)) c("| Gene | Direction | Selected in folds | Frequency | In candidate list (frequency) |", "|---|---|---|---|---|",
      sprintf("| %s | %s | %d of %d | %.2f | %.2f |", fr$gene_id, fr$direction, as.integer(fr$n_selected), nrun, fr$frequency, fr$candidate_frequency)[seq_len(min(nrow(fr), 20))]) else "No gene was selected.", "",
    "## Method details", "",
    "DESeq2 settings are those of conf/thresholds.config (de_padj, de_min_count, de_min_samples; 'auto' = smallest group among the TRAINING samples of the fold). Seeds are fixed: outer folds use seed + 1000 x permutation + repeat, labels of permutation b are drawn with seed + 7919 x b, the model uses the pipeline's lasso_seed and cutoff_seed. Results are deterministic.",
    sprintf("Versions: R %s, glmnet %s, pROC %s.", getRversion(), packageVersion("glmnet"), packageVersion("pROC")), "",
    "## What this does not show", "", LIMIT)
  writeLines(md[!vapply(md, is.null, TRUE)], f_md)
  cat("nested validation written to", normalizePath(OUT), "\n") }
# ---------------------------------------------------------------- dispatch
if (STAGE == "screen") { stage_screen(); quit(save = "no", status = 0) }
if (STAGE == "model") { r <- tryCatch({ stage_model(); NULL }, error = function(e) e); if (inherits(r, "error")) not_run(paste("error in the model stage:", conditionMessage(r))); quit(save = "no", status = 0) }
# ---- STAGE all (orchestrator; base R only)
chk <- tryCatch({
  S <- load_inputs()
  if (S$N < 10) stop(sprintf("only %d samples (at least 10 needed)", S$N))
  if (min(table(S$y)) < 5) stop(sprintf("a class has fewer than 5 samples (%s %d, %s %d)", S$base, sum(S$y == 0), S$pos, sum(S$y == 1)))
  pl <- build_plan(S)
  need <- c(DGE = RS_DGE, CLF = RS_CLF); miss <- need[!file.exists(need)]; if (length(miss)) stop("Rscript not found: ", paste(miss, collapse = ", "))
  pk <- function(rs, pkgs) { o <- suppressWarnings(system2(rs, c("-e", shQuote(sprintf("cat(all(vapply(c(%s), requireNamespace, TRUE, quietly = TRUE)))", paste0("'", pkgs, "'", collapse = ",")))), stdout = TRUE, stderr = FALSE)); isTRUE(any(grepl("TRUE$", o))) }
  if (!pk(RS_DGE, c("DESeq2", "apeglm"))) stop("DESeq2 and apeglm are not both available in ", RS_DGE)
  if (!pk(RS_CLF, c("glmnet", "pROC"))) stop("glmnet and pROC are not both available in ", RS_CLF)
  list(S = S, pl = pl) }, error = function(e) e)
if (inherits(chk, "error")) { if (CHECK) { cat("CHECK FAILED:", conditionMessage(chk), "\n"); quit(save = "no", status = 1) }; not_run(conditionMessage(chk)) }
nt <- length(chk$pl$tasks); nsk <- sum(vapply(chk$pl$tasks, function(t) t$skip != "", TRUE))
if (CORES_REQ > MAXC) cat(sprintf("NOTE: --cores %d reduced to %d (NESTED_MAX_CORES; shared pod)\n", CORES_REQ, MAXC))
cat(sprintf("plan: %d samples, %d outer folds x %d repeats, %d permutations x %d repeats; %d tasks, %d skipped by the plan, %d DESeq2 fits; cores %d%s\n", chk$S$N, KO, NREP, NPERM, PREP, nt, nsk, if (OUTSIDE) 1L else nt - nsk, CORES, if (OUTSIDE) " (comparator: screen once on all samples)" else ""))
if (CHECK) { cat("CHECK OK: inputs readable, both Rscript environments have their packages\n"); quit(save = "no", status = 0) }
passthru <- unlist(lapply(setdiff(names(opt), c("stage", "cores")), function(n) c(paste0("--", n), opt[[n]])))
env1 <- c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1")
run_stage <- function(rs, st) { t0 <- Sys.time(); st_ <- system2(rs, c(shQuote(normalizePath(SELF)), passthru, "--stage", st, "--cores", CORES), env = env1)
  cat(sprintf("stage %s finished, exit %d, %.0f s\n", st, st_, as.numeric(difftime(Sys.time(), t0, units = "secs")))); st_ }
if (!file.exists(file.path(WORK, ".screen_complete"))) {
  dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
  if (run_stage(RS_DGE, "screen") != 0) not_run("the DESeq2 screen stage (r_dge) exited with an error")
  writeLines("done", file.path(WORK, ".screen_complete")) }
if (run_stage(RS_CLF, "model") != 0) not_run("the model stage (r_classifier) exited with an error")
if (!KEEP) unlink(WORK, recursive = TRUE)
quit(save = "no", status = 0)
