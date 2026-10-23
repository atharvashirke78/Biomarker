#!/usr/bin/env Rscript
# Apply the LOCKED diagnostic model (tag lock-diagnostic-v1) to a featureCounts count matrix. Nothing is fitted or tuned here.
#   score(sample) = mean(log2(CPM + 1) over the UP genes) - mean(log2(CPM + 1) over the DOWN genes)
#   CPM = featureCounts count / the sample's featureCounts assigned total (column sum of the matrix) * 1e6
#   call tumour if score > cutoff (strictly greater)        -- lock/diagnostic_v1/final_model_spec.txt
# Genes, signs and cutoff are READ from lock/diagnostic_v1/ after their sha256 is checked against sha256sums_final.txt.
#
# usage (project root, r_classifier env):
#   env/prefix/r_classifier/bin/Rscript scripts/apply_locked_model.R --counts <counts.tsv> --out-dir <new dir>
#       --cohort development|validation [--sample-sheet <sheet.tsv>] [--qc-table <qc.tsv>] [--evaluate]
#   --counts        genes x samples, first column gene_id (versioned Ensembl IDs as in results/04_counts/counts_featurecounts.tsv
#                   or the QUANT contract counts_matrix.tsv)
#   --sample-sheet  contract sheet (sample_id, group, baseline, study, ...); needed for --evaluate
#   --qc-table      sample QC table (srr or sample_id + status); samples with status FAIL are excluded automatically
#                   (locked QC exclusion rule); WARN samples are kept and reported
#   --cohort        label written into every output; 'development' results are always labelled optimistic, not validation
#   --evaluate      compute calls vs labels: errors per class (and per study), the locked pass rule, AUC with DeLong CI.
#                   For --cohort validation this needs --qc-table and a sheet that covers every scored sample.
# Outputs (the out dir must not exist): scores.tsv, excluded_samples.tsv, run_info.tsv [, evaluation.tsv, evaluation_sens_spec.tsv]
suppressPackageStartupMessages(library(pROC))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(evaluate = FALSE)
i <- 1
while (i <= length(args)) {
  k <- sub("^--", "", args[i])
  if (k == "evaluate") { opt$evaluate <- TRUE; i <- i + 1; next }
  if (i + 1 > length(args)) stop("missing value for --", k)
  opt[[gsub("-", "_", k)]] <- args[i + 1]; i <- i + 2
}
die <- function(...) { message("ERROR: ", ...); quit(status = 2) }
if (is.null(opt$counts) || is.null(opt$out_dir) || is.null(opt$cohort)) die("--counts, --out-dir and --cohort are required")
if (!opt$cohort %in% c("development", "validation")) die("--cohort must be development or validation")
if (file.exists(opt$out_dir)) die("output folder exists (never overwritten): ", opt$out_dir)
LOCK <- "lock/diagnostic_v1"

# ---- 1. lock integrity, then read genes / signs / cutoff from the locked files
sha <- function(p) sub(" .*", "", system2("sha256sum", p, stdout = TRUE))
ref <- read.table(file.path(LOCK, "sha256sums_final.txt"), col.names = c("sha", "path"), stringsAsFactors = FALSE)
for (f in c("final_model_genes.tsv", "final_model_spec.txt")) {
  p <- file.path(LOCK, f); want <- ref$sha[ref$path == p]
  if (length(want) != 1 || sha(p) != want) die("lock file does not match sha256sums_final.txt: ", p)
}
genes <- read.delim(file.path(LOCK, "final_model_genes.tsv"), stringsAsFactors = FALSE)
genes$sign <- as.numeric(genes$sign)
spec <- readLines(file.path(LOCK, "final_model_spec.txt"))
cl <- grep("^cutoff = ", spec, value = TRUE)
if (length(cl) != 1) die("cannot read a single 'cutoff = ' line from final_model_spec.txt")
CUT <- as.numeric(sub("^cutoff = ([0-9.eE+-]+).*", "\\1", cl))
up <- genes$gene_id[genes$sign > 0]; dn <- genes$gene_id[genes$sign < 0]
if (!is.finite(CUT) || length(up) < 1 || length(dn) < 1 || length(up) + length(dn) != nrow(genes)) die("unexpected locked model content")

# ---- 2. counts, QC exclusion
cn <- as.matrix(read.delim(opt$counts, row.names = 1, check.names = FALSE))
if (any(duplicated(rownames(cn)))) die("duplicate gene IDs in the counts matrix")
miss <- setdiff(genes$gene_id, rownames(cn))
if (length(miss)) die("locked model genes missing from the counts matrix (versioned IDs must match): ", paste(miss, collapse = ", "))
excl <- data.frame(sample = character(0), reason = character(0))
warn <- character(0)
if (!is.null(opt$qc_table)) {
  qc <- read.delim(opt$qc_table, stringsAsFactors = FALSE)
  idc <- intersect(c("sample_id", "srr"), names(qc))[1]
  if (is.na(idc) || !"status" %in% names(qc)) die("qc table needs sample_id (or srr) and status columns")
  rs <- if ("reason" %in% names(qc)) qc$reason else rep("", nrow(qc))
  fail <- qc[[idc]][qc$status == "FAIL"]
  warn <- intersect(qc[[idc]][qc$status == "WARN"], colnames(cn))
  drop <- intersect(fail, colnames(cn))
  if (length(fail)) excl <- data.frame(sample = fail, reason = paste0("QC FAIL", ifelse(fail %in% drop, " (dropped here)", " (already absent from the counts matrix)"),
                                                                       ": ", rs[match(fail, qc[[idc]])]))
  cn <- cn[, setdiff(colnames(cn), drop), drop = FALSE]
  nq <- setdiff(colnames(cn), qc[[idc]])
  if (length(nq)) die("samples without a QC row: ", paste(nq, collapse = ", "))
}
if (ncol(cn) < 1) die("no samples left to score")

# ---- 3. score (locked formula)
tot <- colSums(cn)
if (any(tot <= 0)) die("a sample has a zero featureCounts total")
L <- log2(sweep(cn[c(up, dn), , drop = FALSE], 2, tot, "/") * 1e6 + 1)
score <- colMeans(L[up, , drop = FALSE]) - colMeans(L[dn, , drop = FALSE])
out <- data.frame(sample = colnames(cn), score = unname(score), call_tumour = as.integer(score > CUT),
                  qc_warn = colnames(cn) %in% warn, cohort = opt$cohort, stringsAsFactors = FALSE)
label <- if (opt$cohort == "development") "development data: optimistic, NOT validation" else "validation data (locked model, applied once)"

# ---- 4. optional evaluation against labels
ev <- NULL
if (opt$evaluate) {
  if (is.null(opt$sample_sheet)) die("--evaluate needs --sample-sheet")
  if (opt$cohort == "validation" && is.null(opt$qc_table)) die("--evaluate on validation data needs --qc-table (locked QC exclusion rule)")
  ss <- read.delim(opt$sample_sheet, stringsAsFactors = FALSE)
  if (!"sample_id" %in% names(ss) || !"group" %in% names(ss)) die("sample sheet needs sample_id and group columns")
  nos <- setdiff(out$sample, ss$sample_id)
  if (length(nos)) die("scored samples missing from the sample sheet: ", paste(nos, collapse = ", "))
  ss <- ss[match(out$sample, ss$sample_id), ]
  pos <- if ("baseline" %in% names(ss)) toupper(ss$baseline) == "FALSE" else ss$group != "normal"
  if (length(unique(ss$group[pos])) != 1 || length(unique(ss$group[!pos])) != 1) die("need exactly one positive (non-baseline) and one baseline group")
  out$truth <- ss$group; out$truth_positive <- as.integer(pos); out$study <- ss$study
  y <- out$truth_positive; p <- out$call_tumour
  err <- function(k) c(errors_positive = sum(y[k] == 1 & p[k] == 0), n_positive = sum(y[k] == 1),
                       errors_baseline = sum(y[k] == 0 & p[k] == 1), n_baseline = sum(y[k] == 0))
  rows <- list(c(scope = "all", err(rep(TRUE, length(y)))))
  st <- unique(out$study)
  for (s in st) rows[[length(rows) + 1]] <- c(scope = paste("study", s), err(out$study == s))
  ev <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  for (cc in names(ev)[-1]) ev[[cc]] <- as.integer(ev[[cc]])
  # locked pass rule: at most 1 error per class, and the same within each study if both studies are present
  scopes <- if (length(st) > 1) seq_len(nrow(ev)) else 1
  ev$within_pass_rule <- ev$errors_positive <= 1 & ev$errors_baseline <= 1
  pass <- all(ev$within_pass_rule[scopes])
  auc_txt <- "n/a (needs both classes)"
  if (length(unique(y)) == 2) {
    r <- roc(y, out$score, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- suppressWarnings(ci.auc(r, method = "delong"))
    auc_txt <- sprintf("%.4f (DeLong 95%% CI %.4f-%.4f)", as.numeric(auc(r)), ci[1], ci[3])
  }
  ev$auc_not_pass_test <- c(auc_txt, rep("", nrow(ev) - 1))
  res <- if (pass) "PASS" else "FAIL"
  if (opt$cohort == "development") res <- paste(res, "(development data only; NOT the validation test)")
  ev$locked_pass_rule_result <- c(res, rep("", nrow(ev) - 1))
  ev$label <- label
}

# ---- 5. write
dir.create(opt$out_dir, recursive = TRUE)
w <- function(x, f) write.table(x, file.path(opt$out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(out, "scores.tsv"); w(excl, "excluded_samples.tsv"); if (!is.null(ev)) w(ev, "evaluation.tsv")
if (!is.null(ev)) {      # NEW file: sensitivity and specificity with Clopper-Pearson 95% intervals (evaluation.tsv is unchanged)
  cp95 <- function(x, n) if (n > 0) { ci <- binom.test(x, n)$conf.int; c(x / n, ci[1], ci[2]) } else c(NA, NA, NA)
  w(do.call(rbind, lapply(seq_len(nrow(ev)), function(i) {
    se <- cp95(ev$n_positive[i] - ev$errors_positive[i], ev$n_positive[i]); sp <- cp95(ev$n_baseline[i] - ev$errors_baseline[i], ev$n_baseline[i])
    data.frame(scope = ev$scope[i], n_positive = ev$n_positive[i], sensitivity = se[1], sens_ci_low = se[2], sens_ci_high = se[3],
               n_baseline = ev$n_baseline[i], specificity = sp[1], spec_ci_low = sp[2], spec_ci_high = sp[3], ci_method = "Clopper-Pearson 95%", label = label) })), "evaluation_sens_spec.tsv")
}
info <- data.frame(key = c("label", "cohort", "lock_tag", "final_model_genes.tsv sha256", "final_model_spec.txt sha256", "cutoff",
                           "n_up", "n_down", "counts", "counts sha256", "sample_sheet", "qc_table", "n_scored", "n_excluded_qc_fail",
                           "n_qc_warn_kept", "evaluate", "R", "pROC", "command"),
                   value = c(label, opt$cohort, "lock-diagnostic-v1", sha(file.path(LOCK, "final_model_genes.tsv")), sha(file.path(LOCK, "final_model_spec.txt")),
                             format(CUT, digits = 12), length(up), length(dn), opt$counts, sha(opt$counts), ifelse(is.null(opt$sample_sheet), "", opt$sample_sheet),
                             ifelse(is.null(opt$qc_table), "", opt$qc_table), nrow(out), nrow(excl), length(warn), opt$evaluate,
                             R.version.string, as.character(packageVersion("pROC")), paste(commandArgs(), collapse = " ")))
w(info, "run_info.tsv")
message(sprintf("scored %d samples (%s); excluded (QC FAIL) %d; cutoff %s; outputs in %s", nrow(out), label, nrow(excl), format(CUT, digits = 12), opt$out_dir))
if (!is.null(ev)) message("locked pass rule: ", ev$locked_pass_rule_result[1], " | AUC (not the pass test): ", ev$auc_not_pass_test[1])
