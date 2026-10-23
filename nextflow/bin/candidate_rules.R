#!/usr/bin/env Rscript
# Candidate rules 1-3b for -entry DISCOVER (Phase 1, 2026-10-15). Generalised COPY of scripts/candidate_selection.R (the original is unchanged).
# Same calls and the same output files; what changed is only where the cohort facts come from:
#   samples, groups, baseline group, batch and patients from the DISCOVER sample sheet (not ncol == 19, 'normal', 'tumour', study IDs);
#   the leave-one-out pass threshold = ceil(robustness_min_fraction x number of refits) (config/questions.yaml); for 19 samples with both routes
#   this must equal the locked loo_min_pass (18, config/candidate_rules.yaml) - checked below, the run stops otherwise.
# Single-route input (no --de-secondary): rule 2 applied to the single route, rule 3a SKIPPED; results are NOT equivalent to the 19-sample method.
# usage: candidate_rules.R --counts <counts.tsv> --sheet <sample_sheet.discover.tsv> --de-primary <deseq2_results_all.tsv> [--de-secondary <deseq2_salmon_results_all.tsv>]
#        --rules config/candidate_rules.yaml --thresholds conf/thresholds.config --questions config/questions.yaml
suppressPackageStartupMessages({ library(DESeq2); library(apeglm) })
a <- commandArgs(trailingOnly = TRUE); opt <- list()
for (i in seq(1, length(a), by = 2)) opt[[sub("^--", "", a[i])]] <- a[i + 1]
for (k in c("counts", "sheet", "de-primary", "rules", "thresholds", "questions")) if (is.null(opt[[k]])) stop("missing --", k)
dual <- !is.null(opt[["de-secondary"]])
# --patient-level-cv true|false (default false = the locked sample-level leave-one-out, even when a patient_id is repeated)
PLC <- tolower(if (is.null(opt[["patient-level-cv"]])) "false" else opt[["patient-level-cv"]]); stopifnot(PLC %in% c("true", "false"))
# --plan-only: TEST HOOK, never passed by the workflow. Needs BIOMARKER_TEST_PLAN_ONLY=1 as well; stops after the refit plan (robustness_refits.tsv + settings line).
PLAN_ONLY <- !is.null(opt[["plan-only"]])
if (PLAN_ONLY && Sys.getenv("BIOMARKER_TEST_PLAN_ONLY") != "1") stop("--plan-only is a test hook and needs BIOMARKER_TEST_PLAN_ONLY=1")
yaml_vals <- function(file, block, key) { l <- readLines(file); i <- grep(paste0("^", block, ":"), l); stopifnot(length(i) == 1)
  j <- i + 1; while (j <= length(l) && grepl("^(\\s|#|$)", l[j])) j <- j + 1; blk <- l[(i + 1):(j - 1)]
  v <- grep(paste0("^\\s*", key, ":"), blk, value = TRUE); as.numeric(sub("^\\s*[a-z0-9_]+:\\s*([0-9.]+).*$", "\\1", v)) }
read_key <- function(file, key, sep) { l <- grep(paste0("^\\s*", key, "\\s*", sep), readLines(file), value = TRUE); stopifnot(length(l) == 1)
  trimws(gsub("['\"]", "", sub("(#|//).*$", "", sub(paste0("^\\s*", key, "\\s*", sep), "", l)))) }
RULES <- opt$rules
P_CUT <- yaml_vals(RULES, "rule1_significance", "padj_max"); LFC_CUT <- yaml_vals(RULES, "rule1_significance", "abs_shrunk_log2fc_min")
bm <- yaml_vals(RULES, "rule2_measurable", "basemean_min"); stopifnot(length(bm) == 2); BM_MIN <- bm[1]; BM_STRICT <- bm[2]
stopifnot(length(P_CUT) == 1, length(LFC_CUT) == 1, is.finite(c(P_CUT, LFC_CUT, BM_MIN, BM_STRICT)))
ALPHA <- as.numeric(read_key(opt$thresholds, "de_padj", "=")); MIN_COUNT <- as.numeric(read_key(opt$thresholds, "de_min_count", "="))
DE_MIN_SAMPLES <- read_key(opt$thresholds, "de_min_samples", "=")
# values resolved by Nextflow (command-line / profile overrides included) win over the text of thresholds.config, so the refits use the same DESeq2 settings as the main fit
if (!is.null(opt[["de-padj"]])) ALPHA <- as.numeric(opt[["de-padj"]])
if (!is.null(opt[["de-min-count"]])) MIN_COUNT <- as.numeric(opt[["de-min-count"]])
if (!is.null(opt[["de-min-samples"]])) DE_MIN_SAMPLES <- opt[["de-min-samples"]]
LOCKED_LOO_MIN_PASS <- as.integer(read_key(RULES, "loo_min_pass", ":"))
METHOD <- read_key(opt$questions, "robustness_method", ":"); MIN_FRAC <- as.numeric(read_key(opt$questions, "robustness_min_fraction", ":"))
N_SUB <- as.integer(read_key(opt$questions, "robustness_subsamples", ":")); SUB_FRAC <- as.numeric(read_key(opt$questions, "robustness_subsample_fraction", ":"))
ROB_SEED <- as.integer(read_key(opt$questions, "robustness_seed", ":"))
stopifnot(is.finite(ALPHA), is.finite(MIN_COUNT), METHOD %in% c("loo", "subsample"))
# ---- OPT-IN settings for rule 3a (config/questions.yaml; every key optional, the defaults reproduce the locked behaviour exactly)
#   rule3a_mode: strict (default) = rule 1 on BOTH routes, same direction (locked rule) | concordance = featureCounts passes rule 1, Salmon has the same direction and
#   padj < rule3a_salmon_padj_max and |shrunk log2FC| >= rule3a_salmon_min_abs_lfc. The two Salmon criteria are PROVISIONAL (no source).
#   rule3a_borderline_padj_margin / rule3a_borderline_lfc_margin: width of the "salmon_borderline" flag around the rule-1 cutoffs (flag only, never changes a call).
opt_key <- function(key, default) { l <- grep(paste0("^\\s*", key, "\\s*:"), readLines(opt$questions), value = TRUE)
  if (length(l) == 0) return(default); if (length(l) > 1) stop("config key ", key, " appears more than once in ", opt$questions)
  v <- trimws(gsub("['\"]", "", sub("#.*$", "", sub(paste0("^\\s*", key, "\\s*:"), "", l)))); if (is.numeric(default)) as.numeric(v) else v }
R3A_MODE <- opt_key("rule3a_mode", "strict"); if (!R3A_MODE %in% c("strict", "concordance")) stop("rule3a_mode must be strict or concordance, not '", R3A_MODE, "'")
SAL_PADJ_MAX <- opt_key("rule3a_salmon_padj_max", 0.1); SAL_MIN_LFC <- opt_key("rule3a_salmon_min_abs_lfc", 0.5)
BORD_P <- opt_key("rule3a_borderline_padj_margin", 0.01); BORD_L <- opt_key("rule3a_borderline_lfc_margin", 0.1)
stopifnot(is.finite(c(SAL_PADJ_MAX, SAL_MIN_LFC, BORD_P, BORD_L)))
if (R3A_MODE != "strict") cat(sprintf("WARNING: rule3a_mode = %s (NOT the locked rule 3a). Salmon route needs same direction, padj < %g and |shrunk log2FC| >= %g (PROVISIONAL, no source). Results are not equivalent to the locked 19-sample method.\n", R3A_MODE, SAL_PADJ_MAX, SAL_MIN_LFC))

# ---- inputs (cohort facts from the sample sheet)
cn <- as.matrix(read.delim(opt$counts, row.names = 1, check.names = FALSE))
sm <- read.delim(opt$sheet, stringsAsFactors = FALSE, colClasses = "character"); rownames(sm) <- sm$sample_id
stopifnot(setequal(colnames(cn), sm$sample_id))
base_group <- unique(sm$group[toupper(sm$baseline) == "TRUE"]); stopifnot(length(base_group) == 1)
ss <- sm[colnames(cn), ]; ss$batch <- factor(ss$batch); ss$group <- relevel(factor(ss$group), ref = base_group)
test_group <- setdiff(levels(ss$group), base_group)
N <- ncol(cn)
design <- if (nlevels(ss$batch) > 1) ~ batch + group else ~ group
print(table(ss$batch, ss$group))
MIN_SAMPLES <- if (DE_MIN_SAMPLES == "auto") min(table(ss$group)) else as.numeric(DE_MIN_SAMPLES)   # from the full cohort, not per refit
# refits: leave-one-out over samples (locked method) or patient-level subsamples
PATIENT_LEVEL <- FALSE
if (METHOD == "loo") {
  # locked method = leave one SAMPLE out. Leave one PATIENT out (all samples of that patient together) only when --patient-level-cv true AND a patient_id is repeated.
  pp <- if ("patient_id" %in% names(ss)) ss[colnames(cn), "patient_id"] else colnames(cn)
  PATIENT_LEVEL <- PLC == "true" && any(duplicated(pp))
  if (!PATIENT_LEVEL && any(duplicated(pp))) message("Repeated patient_id found; folds and refits are sample-level (the locked method). Use --patient_level_cv for patient-level.")
  refits <- if (PATIENT_LEVEL) setNames(lapply(unique(pp), function(p) setdiff(colnames(cn), colnames(cn)[pp == p])), unique(pp))
            else setNames(lapply(colnames(cn), function(s) setdiff(colnames(cn), s)), colnames(cn))
} else {
  set.seed(ROB_SEED); pats <- unique(ss$patient_id)
  refits <- setNames(lapply(seq_len(N_SUB), function(r) { keep <- sample(pats, ceiling(SUB_FRAC * length(pats))); colnames(cn)[ss$patient_id %in% keep] }), paste0("subsample_", seq_len(N_SUB)))
}
drop_txt <- function(s) if (METHOD != "loo") paste(s, length(refits_all[[s]]), "samples") else if (PATIENT_LEVEL) paste0("dropping patient ", s, " (samples ", paste(setdiff(colnames(cn), refits_all[[s]]), collapse = ", "), ")") else paste0("dropping ", s, " (", ss[s, "label"], ")")
refits_all <- refits
# A refit is SKIPPED (not stopped) when dropping its sample(s) would leave a batch with only one group: the batch + group design is then not estimable.
# Same condition as the stopifnot inside fit(), which stays as a safety net. With the 19 development samples no refit is ever skipped.
MIN_REFITS <- 10L        # PROVISIONAL: below this the "90 % of refits" rule cannot tolerate one failed refit (ceil(0.9 x 9) = 9 = all); optional key robustness_min_refits in questions.yaml
mr <- grep("^\\s*robustness_min_refits\\s*:", readLines(opt$questions), value = TRUE); if (length(mr) == 1) MIN_REFITS <- as.integer(sub("^[^:]*:\\s*([0-9]+).*$", "\\1", mr))
refit_status <- do.call(rbind, lapply(names(refits), function(s) {
  d0 <- ss[refits[[s]], , drop = FALSE]; tb <- table(droplevels(d0$batch), droplevels(d0$group))
  miss <- which(tb == 0, arr.ind = TRUE)
  note <- if (nrow(miss)) paste0("skipped: batch ", rownames(tb)[miss[1, 1]], " would have only one group (no ", colnames(tb)[miss[1, 2]], " sample left)") else ""
  data.frame(refit = s, n_samples = length(refits[[s]]), status = if (nrow(miss)) "skipped" else "run", note = note, stringsAsFactors = FALSE) }))
N_REFIT_ALL <- nrow(refit_status); N_SKIPPED <- sum(refit_status$status == "skipped")
# H2d: feasibility of rule 3b (written before any stop so the file exists even when no refit can run)
write.table(data.frame(n_samples = N, n_refits_planned = N_REFIT_ALL, n_refits_run = N_REFIT_ALL - N_SKIPPED, n_skipped = N_SKIPPED, min_refits_required = MIN_REFITS,
  evaluable = (N_REFIT_ALL - N_SKIPPED) >= MIN_REFITS), "robustness_feasibility.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
for (i in which(refit_status$status == "skipped")) message("robustness refit SKIPPED: ", drop_txt(refit_status$refit[i]), ": ", sub("^skipped: ", "", refit_status$note[i]))
write.table(refit_status, "robustness_refits.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
refits <- refits[refit_status$status == "run"]
N_REFIT <- length(refits)
if (N_REFIT == 0) stop("all ", N_REFIT_ALL, " robustness refits were skipped (every one would leave a batch with only one group): rule 3b cannot be evaluated. See robustness_refits.tsv")
ROB_EVALUABLE <- N_REFIT >= MIN_REFITS
LOO_MIN_PASS <- as.integer(ceiling(MIN_FRAC * N_REFIT - 1e-9))     # scaled to the refits that ran; with all 19 refits: 18, the locked value
if (METHOD == "loo" && !PATIENT_LEVEL && N == 19 && dual && LOO_MIN_PASS != LOCKED_LOO_MIN_PASS) stop("derived loo_min_pass ", LOO_MIN_PASS, " != locked ", LOCKED_LOO_MIN_PASS, " for the 19-sample method")
if (N_SKIPPED > 0) cat(sprintf("robustness: %d of %d refits run, %d skipped (see robustness_refits.tsv)%s\n", N_REFIT, N_REFIT_ALL, N_SKIPPED, if (ROB_EVALUABLE) "" else sprintf("; fewer than %d refits ran: rule 3b NOT EVALUABLE", MIN_REFITS)))
cat(sprintf("settings: alpha %g | pre-filter >= %g counts in >= %g samples (de_min_samples = %s, from full cohort) | robustness %s, %d refits, pass >= %d (= ceil(%g x %d)) | routes: %s\n",
            ALPHA, MIN_COUNT, MIN_SAMPLES, DE_MIN_SAMPLES, METHOD, N_REFIT, LOO_MIN_PASS, MIN_FRAC, N_REFIT, if (dual) "featureCounts + Salmon" else "single route"))
if (PLAN_ONLY) quit(save = "no", status = 0)
fcres <- read.delim(opt[["de-primary"]], stringsAsFactors = FALSE)
sares <- if (dual) read.delim(opt[["de-secondary"]], stringsAsFactors = FALSE) else NULL

# ---- one fit = the cohort call (deseq2.R design), applied to a subset of columns
fit <- function(cols) {
  d <- DESeqDataSetFromMatrix(round(cn[, cols, drop = FALSE]), colData = ss[cols, ], design = design)
  d$batch <- droplevels(d$batch); d$group <- droplevels(d$group)
  stopifnot(all(table(d$batch, d$group) > 0))
  d <- d[rowSums(counts(d) >= MIN_COUNT) >= MIN_SAMPLES, ]
  d <- DESeq(d, quiet = TRUE)
  coef <- paste0("group_", test_group, "_vs_", base_group); stopifnot(coef %in% resultsNames(d))
  res <- results(d, name = coef, alpha = ALPHA); shr <- lfcShrink(d, coef = coef, type = "apeglm", res = res, quiet = TRUE)
  data.frame(gene_id = rownames(res), baseMean = res$baseMean, lfc = res$log2FoldChange, lfc_shr = shr$log2FoldChange, padj = res$padj, stringsAsFactors = FALSE)
}
r1 <- function(padj, l) !is.na(padj) & padj < P_CUT & !is.na(l) & abs(l) >= LFC_CUT

# ---- full-cohort refit (must reproduce the primary DE table), then the robustness refits
full <- fit(colnames(cn)); rownames(full) <- full$gene_id
m <- merge(full, fcres, by = "gene_id", suffixes = c("", "_old"))
chk <- data.frame(n_genes_new = nrow(full), n_genes_old = nrow(fcres), n_compared = nrow(m),
  max_abs_diff_log2FC = max(abs(m$lfc - m$log2FoldChange), na.rm = TRUE), max_abs_diff_shrunk = max(abs(m$lfc_shr - m$log2FoldChange_shrunk), na.rm = TRUE),
  max_abs_diff_padj = max(abs(m$padj - m$padj_old), na.rm = TRUE), n_rule1_new = sum(r1(full$padj, full$lfc_shr)), n_rule1_old = sum(r1(fcres$padj, fcres$log2FoldChange_shrunk)))
print(chk); write.table(chk, "reproduction_check_19_sample_fit.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
full_dir <- sign(full$lfc_shr); names(full_dir) <- full$gene_id
genes <- full$gene_id
loo_pass <- matrix(FALSE, nrow = length(genes), ncol = N_REFIT, dimnames = list(genes, names(refits)))
for (s in names(refits)) {
  message("robustness refit: ", drop_txt(s))
  f <- fit(refits[[s]]); ok <- r1(f$padj, f$lfc_shr) & sign(f$lfc_shr) == full_dir[f$gene_id]
  loo_pass[f$gene_id[ok], s] <- TRUE
}
write.table(data.frame(gene_id = genes, loo_pass, check.names = FALSE), "loo_pass_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- candidate table: every gene that passes rule 1 on either route, one row per gene_id
fcres$r1 <- r1(fcres$padj, fcres$log2FoldChange_shrunk)
if (dual) sares$r1 <- r1(sares$padj, sares$log2FoldChange_shrunk)
ids <- if (dual) union(fcres$gene_id[fcres$r1], sares$gene_id[sares$r1]) else fcres$gene_id[fcres$r1]
f <- fcres[match(ids, fcres$gene_id), ]
s <- if (dual) sares[match(ids, sares$gene_id), ] else data.frame(gene_name = rep(NA, length(ids)), gene_type = NA, baseMean = NA, log2FoldChange_shrunk = NA, log2FoldChange = NA, padj = NA, r1 = NA)
t <- data.frame(gene_id = ids, symbol = ifelse(!is.na(f$gene_name), f$gene_name, s$gene_name), gene_biotype = ifelse(!is.na(f$gene_type), f$gene_type, s$gene_type),
  baseMean_fc = f$baseMean, baseMean_salmon = s$baseMean, log2FC_shrunk_fc = f$log2FoldChange_shrunk, log2FC_shrunk_salmon = s$log2FoldChange_shrunk,
  log2FC_unshrunk_fc = f$log2FoldChange, log2FC_unshrunk_salmon = s$log2FoldChange, padj_fc = f$padj, padj_salmon = s$padj, stringsAsFactors = FALSE)
t$rule1_fc <- ifelse(is.na(f$r1), FALSE, f$r1); t$rule1_salmon <- ifelse(is.na(s$r1), FALSE, s$r1)
t$direction_fc <- ifelse(is.na(t$log2FC_shrunk_fc), NA, ifelse(t$log2FC_shrunk_fc > 0, "up", "down"))
t$direction_salmon <- ifelse(is.na(t$log2FC_shrunk_salmon), NA, ifelse(t$log2FC_shrunk_salmon > 0, "up", "down"))
t$rule1_pass <- t$rule1_fc | t$rule1_salmon
if (dual) {
  t$rule2_pass <- !is.na(t$baseMean_fc) & !is.na(t$baseMean_salmon) & t$baseMean_fc >= BM_MIN & t$baseMean_salmon >= BM_MIN
  t$baseMean_ge_50_both <- !is.na(t$baseMean_fc) & !is.na(t$baseMean_salmon) & t$baseMean_fc >= BM_STRICT & t$baseMean_salmon >= BM_STRICT
  same_dir <- !is.na(t$direction_fc) & !is.na(t$direction_salmon) & t$direction_fc == t$direction_salmon
  t$rule3a_pass <- if (R3A_MODE == "strict") t$rule1_fc & t$rule1_salmon & same_dir      # locked rule
    else t$rule1_fc & same_dir & !is.na(t$padj_salmon) & t$padj_salmon < SAL_PADJ_MAX & !is.na(t$log2FC_shrunk_salmon) & abs(t$log2FC_shrunk_salmon) >= SAL_MIN_LFC
} else {   # single route: rule 2 on the one route (MODIFIED), rule 3a SKIPPED (counted as passed for all_rules_pass, never reported as passed)
  t$rule2_pass <- !is.na(t$baseMean_fc) & t$baseMean_fc >= BM_MIN
  t$baseMean_ge_50_both <- !is.na(t$baseMean_fc) & t$baseMean_fc >= BM_STRICT
  t$rule3a_pass <- NA
}
lp <- loo_pass[match(ids, rownames(loo_pass)), , drop = FALSE]; tested <- !is.na(match(ids, rownames(loo_pass)))
t$loo_n_pass <- ifelse(tested, rowSums(lp), NA); t$loo_fraction <- t$loo_n_pass / N_REFIT; t$loo_n_fail <- N_REFIT - t$loo_n_pass
t$loo_failed_samples <- vapply(seq_along(ids), function(i) if (!tested[i]) "not tested in 19-sample fit" else paste(colnames(lp)[!lp[i, ]], collapse = ";"), "")
t$rule3b_pass <- if (ROB_EVALUABLE) !is.na(t$loo_n_pass) & t$loo_n_pass >= LOO_MIN_PASS else rep(NA, nrow(t))      # NA = NOT EVALUABLE
t$all_rules_pass <- t$rule1_pass & t$rule2_pass & (if (dual) t$rule3a_pass else TRUE) & (if (ROB_EVALUABLE) t$rule3b_pass else FALSE)
why <- function(i) { r <- character(0)
  if (!t$rule2_pass[i]) r <- c(r, if (dual) sprintf("rule2: baseMean fc %s / salmon %s (need >= %g on both)", format(t$baseMean_fc[i], digits = 3), format(t$baseMean_salmon[i], digits = 3), BM_MIN)
                                  else sprintf("rule2 (single route): baseMean %s (need >= %g)", format(t$baseMean_fc[i], digits = 3), BM_MIN))
  if (dual && !t$rule3a_pass[i]) r <- c(r, if (!t$rule1_fc[i]) "rule3a: fails rule 1 on featureCounts"
    else if (R3A_MODE == "strict") { if (!t$rule1_salmon[i]) "rule3a: fails rule 1 on Salmon" else "rule3a: opposite direction" }
    else if (is.na(t$direction_salmon[i])) "rule3a (concordance): no Salmon result"
    else if (t$direction_fc[i] != t$direction_salmon[i]) "rule3a (concordance): opposite direction"
    else sprintf("rule3a (concordance): Salmon padj %s / |shrunk log2FC| %s not within padj < %g and >= %g", format(t$padj_salmon[i], digits = 3), format(abs(t$log2FC_shrunk_salmon[i]), digits = 3), SAL_PADJ_MAX, SAL_MIN_LFC))
  if (!ROB_EVALUABLE) r <- c(r, sprintf("rule3b: NOT EVALUABLE (only %d of %d refits ran, %d skipped; at least %d needed)", N_REFIT, N_REFIT_ALL, N_SKIPPED, MIN_REFITS))
  else if (!t$rule3b_pass[i]) r <- c(r, sprintf("rule3b: passes %s of %d leave-one-out refits (need >= %d)", t$loo_n_pass[i], N_REFIT, LOO_MIN_PASS))
  paste(r, collapse = "; ") }
t$fail_reason <- vapply(seq_len(nrow(t)), why, "")
# NEW columns (always written, after the existing ones; existing columns and values are unchanged)
t$salmon_borderline <- if (dual) !is.na(t$padj_salmon) & !is.na(t$log2FC_shrunk_salmon) & (abs(t$padj_salmon - P_CUT) <= BORD_P + 1e-12 | abs(abs(t$log2FC_shrunk_salmon) - LFC_CUT) <= BORD_L + 1e-12) else rep(NA, nrow(t))
t$rule3a_mode_used <- if (dual) R3A_MODE else "not_applicable_single_route"
t <- t[order(!t$all_rules_pass, -t$loo_n_pass, t$padj_fc), ]
write.table(t, "candidate_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- funnel (sequential), all genes and protein-coding
pc <- t$gene_biotype == "protein_coding"
step <- function(n, k) data.frame(step = n, all_genes = sum(k, na.rm = TRUE), protein_coding = sum(k & pc, na.rm = TRUE))
if (dual) {
  fun <- rbind(step("rule 1 on featureCounts", t$rule1_fc), step("rule 1 on Salmon", t$rule1_salmon), step("rule 1 on either route (candidate table)", t$rule1_pass),
    step(sprintf("+ rule 2 (baseMean >= %g on both routes)", BM_MIN), t$rule1_pass & t$rule2_pass), step(if (R3A_MODE == "strict") "+ rule 3a (rule 1 on both routes, same direction)" else sprintf("+ rule 3a CONCORDANCE (featureCounts rule 1; Salmon same direction, padj < %g, |shrunk log2FC| >= %g)", SAL_PADJ_MAX, SAL_MIN_LFC), t$rule1_pass & t$rule2_pass & t$rule3a_pass),
    step(if (ROB_EVALUABLE) sprintf("+ rule 3b (LOO >= %d of %d, featureCounts)", LOO_MIN_PASS, N_REFIT) else sprintf("+ rule 3b NOT EVALUABLE (%d of %d refits ran; at least %d needed): no candidate declared", N_REFIT, N_REFIT_ALL, MIN_REFITS), t$all_rules_pass), step(sprintf("  of which baseMean >= %g on both (stricter tier)", BM_STRICT), t$all_rules_pass & t$baseMean_ge_50_both))
} else {
  fun <- rbind(step("rule 1 (single route)", t$rule1_fc),
    step(sprintf("+ rule 2 MODIFIED: baseMean >= %g on the single route (defined for both routes)", BM_MIN), t$rule1_pass & t$rule2_pass),
    step("+ rule 3a SKIPPED (needs both routes)", t$rule1_pass & t$rule2_pass),
    step(if (ROB_EVALUABLE) sprintf("+ rule 3b (%s, >= %d of %d refits)", METHOD, LOO_MIN_PASS, N_REFIT) else sprintf("+ rule 3b NOT EVALUABLE (%d of %d refits ran; at least %d needed): no candidate declared", N_REFIT, N_REFIT_ALL, MIN_REFITS), t$all_rules_pass), step(sprintf("  of which baseMean >= %g (stricter tier, single route)", BM_STRICT), t$all_rules_pass & t$baseMean_ge_50_both))
}
write.table(fun, "funnel_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE); print(fun)

# ---- per-dropped-sample failures (leave-one-out only)
if (METHOD == "loo") {
  unit_of <- as.list(setNames(if (PATIENT_LEVEL) ss[colnames(cn), "patient_id"] else colnames(cn), colnames(cn)))   # the refit (column) in which each sample was left out
  fc1 <- fcres$gene_id[fcres$r1]; k <- t$gene_id[t$rule1_pass & t$rule2_pass & (if (dual) t$rule3a_pass else TRUE)]
  bys <- data.frame(srr = colnames(cn), label = ss$label, study = as.character(ss$study), group = as.character(ss$group),
    fail_among_fc_rule1_genes = sapply(colnames(cn), function(x) { u <- unit_of[[x]]; if (u %in% colnames(loo_pass)) sum(!loo_pass[fc1, u]) else NA_integer_ }), n_fc_rule1_genes = length(fc1),
    fail_among_rule123a_genes = sapply(colnames(cn), function(x) { u <- unit_of[[x]]; if (u %in% colnames(loo_pass)) sum(!loo_pass[k, u]) else NA_integer_ }), n_rule123a_genes = length(k), row.names = NULL)
  write.table(bys, "loo_failures_by_dropped_sample.tsv", sep = "\t", quote = FALSE, row.names = FALSE); print(bys)
}
note <- c(if (PATIENT_LEVEL) "robustness: leave-one-PATIENT-out (a patient_id is repeated, so all samples of a patient are left out together)",
  sprintf("routes: %s", if (dual) "featureCounts + Salmon (both rules 2 and 3a as defined in config/candidate_rules.yaml)" else "SINGLE ROUTE"),
  if (!dual) c("rule 2: MODIFIED (applied to the single route; the rule is defined for both routes)", "rule 3a: SKIPPED (needs both quantification routes)",
               "These candidates are NOT equivalent to the 19-sample method (which used both routes)."),
  sprintf("rule 3b: %s, %d refits, pass >= %d = ceil(%g x %d); locked loo_min_pass for the 19-sample method = %d", METHOD, N_REFIT, LOO_MIN_PASS, MIN_FRAC, N_REFIT, LOCKED_LOO_MIN_PASS),
  if (N_SKIPPED > 0) sprintf("robustness refits: %d of %d run, %d SKIPPED (a batch would have only one group); the rule 3b threshold is ceil(%g x %d) on the refits that ran. Details: robustness_refits.tsv", N_REFIT, N_REFIT_ALL, N_SKIPPED, MIN_FRAC, N_REFIT),
  if (!ROB_EVALUABLE) sprintf("RULE 3B NOT EVALUABLE: only %d refits ran (at least %d needed); no candidate can be declared and no model can be fitted from this run", N_REFIT, MIN_REFITS),
  sprintf("rule 3b feasibility: %d samples, %d refits planned, %d run, %d skipped, at least %d required: rule 3b is %s (robustness_feasibility.tsv)", N, N_REFIT_ALL, N_REFIT, N_SKIPPED, MIN_REFITS, if (ROB_EVALUABLE) "evaluable" else "not evaluable"),
  if (!ROB_EVALUABLE) sprintf("rule 3b is not evaluable because only %d of %d planned refits could be run (%d skipped because dropping the sample would leave a batch with one group only; the minimum is %d because below it the %g leave-one-out fraction tolerates no failed refit; PROVISIONAL). With %d samples this is a property of the cohort size and batch layout, not evidence about the genes: more samples per group and batch are needed.", N_REFIT, N_REFIT_ALL, N_SKIPPED, MIN_REFITS, MIN_FRAC, N),
  if (dual) sprintf("rule 3a mode: %s%s", R3A_MODE, if (R3A_MODE == "strict") " (locked rule: rule 1 on both routes, same direction)" else sprintf(" (NOT the locked rule: featureCounts passes rule 1; Salmon same direction, padj < %g, |shrunk log2FC| >= %g; both Salmon criteria PROVISIONAL, no source)", SAL_PADJ_MAX, SAL_MIN_LFC))
  else "rule 3a mode: not applicable (single route)",
  if (R3A_MODE != "strict" && dual) "WARNING: rule3a_mode is not strict; these candidates are not equivalent to the locked 19-sample method",
  if (dual) sprintf("Salmon borderline genes (padj within %g of %g or |shrunk log2FC| within %g of %g): %d of %d genes in the candidate table, %d of %d genes passing rules 1 and 2; the Salmon route is not exactly repeatable (run-to-run differences of Salmon)", BORD_P, P_CUT, BORD_L, LFC_CUT, sum(t$salmon_borderline), nrow(t), sum(t$salmon_borderline & t$rule1_pass & t$rule2_pass), sum(t$rule1_pass & t$rule2_pass)),
  sprintf("candidates passing all applied rules: %d", sum(t$all_rules_pass)))
writeLines(note, "candidate_rules_method_note.txt")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
writeLines(c(paste("R", getRversion()), paste("DESeq2", packageVersion("DESeq2")), paste("apeglm", packageVersion("apeglm"))), "tool_versions.txt")
