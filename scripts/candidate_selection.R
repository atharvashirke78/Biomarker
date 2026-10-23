#!/usr/bin/env Rscript
# Candidate selection, rules 1-3 (config/candidate_rules.yaml). Run from the project root:
#   env/prefix/r_dge/bin/Rscript scripts/candidate_selection.R <counts_dir> <dge_dir> <out_dir> [--force]
# No defaults. counts_dir holds counts_featurecounts.tsv (sample_map.tsv is read from its parent folder); dge_dir holds the 19-sample deseq2_results_all.tsv and deseq2_salmon_results_all.tsv.
# Reads (never writes) those; writes only into out_dir. Refuses to write if out_dir already contains files, unless --force is given.
# Does not read the 15 validation samples.
# The leave-one-out refit repeats the call of nextflow/bin/deseq2.R (same design, pre-filter, results(alpha = de_padj), lfcShrink apeglm).
# rc1-prep: alpha (de_padj), pre-filter count (de_min_count) and n (de_min_samples) are READ from conf/thresholds.config (the file deseq2.R's parameters come from); loo_min_pass from config/candidate_rules.yaml. 'auto' n = smallest group of the FULL 19-sample cohort, fixed for all 19 refits.
suppressPackageStartupMessages({ library(DESeq2); library(apeglm) })
args <- commandArgs(trailingOnly = TRUE); force <- "--force" %in% args; args <- setdiff(args, "--force")
if (length(args) != 3) stop("usage: candidate_selection.R <counts_dir> <dge_dir> <out_dir> [--force]  (no defaults)")
counts_dir <- args[1]; dge_dir <- args[2]; out <- args[3]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) > 0 && !force) stop("ERROR: output folder ", out, " already contains files; nothing written. Use a new folder (or --force).")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
yaml_vals <- function(file, block, key) { l <- readLines(file); i <- grep(paste0("^", block, ":"), l); stopifnot(length(i) == 1)
  j <- i + 1; while (j <= length(l) && grepl("^(\\s|#|$)", l[j])) j <- j + 1; blk <- l[(i + 1):(j - 1)]
  v <- grep(paste0("^\\s*", key, ":"), blk, value = TRUE); as.numeric(sub("^\\s*[a-z0-9_]+:\\s*([0-9.]+).*$", "\\1", v)) }   # numbers of a key inside a top-level yaml block (in file order)
RULES <- "config/candidate_rules.yaml"                                           # rules 1-2 and the baseMean-50 flag are read from here, not typed
P_CUT <- yaml_vals(RULES, "rule1_significance", "padj_max"); LFC_CUT <- yaml_vals(RULES, "rule1_significance", "abs_shrunk_log2fc_min")
bm <- yaml_vals(RULES, "rule2_measurable", "basemean_min"); stopifnot(length(bm) == 2); BM_MIN <- bm[1]; BM_STRICT <- bm[2]   # rule 2, stricter tier (flag only)
stopifnot(length(P_CUT) == 1, length(LFC_CUT) == 1, is.finite(c(P_CUT, LFC_CUT, BM_MIN, BM_STRICT)))
read_key <- function(file, key, sep) { l <- grep(paste0("^\\s*", key, "\\s*", sep), readLines(file), value = TRUE); stopifnot(length(l) == 1)
  trimws(gsub("['\"]", "", sub("(#|//).*$", "", sub(paste0("^\\s*", key, "\\s*", sep), "", l)))) }
ALPHA <- as.numeric(read_key("conf/thresholds.config", "de_padj", "="))                 # same parameter deseq2.R receives
MIN_COUNT <- as.numeric(read_key("conf/thresholds.config", "de_min_count", "="))
DE_MIN_SAMPLES <- read_key("conf/thresholds.config", "de_min_samples", "=")              # 'auto' or a number
LOO_MIN_PASS <- as.integer(read_key("config/candidate_rules.yaml", "loo_min_pass", ":"))
stopifnot(is.finite(ALPHA), is.finite(MIN_COUNT), is.finite(LOO_MIN_PASS))

# ---- inputs
cn <- as.matrix(read.delim(file.path(counts_dir, "counts_featurecounts.tsv"), row.names = 1, check.names = FALSE))
sm <- read.delim(file.path(dirname(normalizePath(counts_dir)), "sample_map.tsv"), stringsAsFactors = FALSE); rownames(sm) <- sm$srr
stopifnot(ncol(cn) == 19, all(colnames(cn) %in% sm$srr))   # deliberate guard: rule 3b (loo_min_pass 18 of 19) is locked for this 19-sample discovery cohort; a new cohort needs a new decision, not a silent re-use
N_LOO <- ncol(cn)
ss <- sm[colnames(cn), ]; ss$study <- factor(ss$study); ss$group <- relevel(factor(ss$group), ref = "normal")
print(table(ss$study, ss$group))
MIN_SAMPLES <- if (DE_MIN_SAMPLES == "auto") min(table(ss$group)) else as.numeric(DE_MIN_SAMPLES)   # from the full cohort, not per refit
cat(sprintf("settings: alpha %g | pre-filter >= %g counts in >= %g samples (de_min_samples = %s, from full cohort) | loo_min_pass %d\n", ALPHA, MIN_COUNT, MIN_SAMPLES, DE_MIN_SAMPLES, LOO_MIN_PASS))
fcres <- read.delim(file.path(dge_dir, "deseq2_results_all.tsv"), stringsAsFactors = FALSE)
sares <- read.delim(file.path(dge_dir, "deseq2_salmon_results_all.tsv"), stringsAsFactors = FALSE)

# ---- one fit = the 19-sample call, applied to a subset of columns
fit <- function(cols) {
  d <- DESeqDataSetFromMatrix(round(cn[, cols, drop = FALSE]), colData = ss[cols, ], design = ~ study + group)
  d$study <- droplevels(d$study); d$group <- droplevels(d$group)
  stopifnot(all(table(d$study, d$group) > 0))
  d <- d[rowSums(counts(d) >= MIN_COUNT) >= MIN_SAMPLES, ]
  d <- DESeq(d, quiet = TRUE)
  coef <- "group_tumour_vs_normal"; stopifnot(coef %in% resultsNames(d))
  res <- results(d, name = coef, alpha = ALPHA); shr <- lfcShrink(d, coef = coef, type = "apeglm", res = res, quiet = TRUE)
  data.frame(gene_id = rownames(res), baseMean = res$baseMean, lfc = res$log2FoldChange, lfc_shr = shr$log2FoldChange, padj = res$padj, stringsAsFactors = FALSE)
}
r1 <- function(padj, l) !is.na(padj) & padj < P_CUT & !is.na(l) & abs(l) >= LFC_CUT

# ---- 19-sample refit (check that this script reproduces the existing featureCounts table), then 19 leave-one-out refits
full <- fit(colnames(cn)); rownames(full) <- full$gene_id
m <- merge(full, fcres, by = "gene_id", suffixes = c("", "_old"))
chk <- data.frame(n_genes_new = nrow(full), n_genes_old = nrow(fcres), n_compared = nrow(m),
  max_abs_diff_log2FC = max(abs(m$lfc - m$log2FoldChange), na.rm = TRUE), max_abs_diff_shrunk = max(abs(m$lfc_shr - m$log2FoldChange_shrunk), na.rm = TRUE),
  max_abs_diff_padj = max(abs(m$padj - m$padj_old), na.rm = TRUE), n_rule1_new = sum(r1(full$padj, full$lfc_shr)), n_rule1_old = sum(r1(fcres$padj, fcres$log2FoldChange_shrunk)))
print(chk); write.table(chk, file.path(out, "reproduction_check_19_sample_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
full_dir <- sign(full$lfc_shr); names(full_dir) <- full$gene_id
genes <- full$gene_id
loo_pass <- matrix(FALSE, nrow = length(genes), ncol = 19, dimnames = list(genes, colnames(cn)))
for (s in colnames(cn)) {
  message("leave-one-out: dropping ", s, " (", ss[s, "label"], ")")
  f <- fit(setdiff(colnames(cn), s)); ok <- r1(f$padj, f$lfc_shr) & sign(f$lfc_shr) == full_dir[f$gene_id]
  loo_pass[f$gene_id[ok], s] <- TRUE
}
write.table(data.frame(gene_id = genes, loo_pass, check.names = FALSE), file.path(out, "loo_pass_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- candidate table: every gene that passes rule 1 on featureCounts or on Salmon (full 19-sample tables), one row per gene_id
fcres$r1 <- r1(fcres$padj, fcres$log2FoldChange_shrunk); sares$r1 <- r1(sares$padj, sares$log2FoldChange_shrunk)
ids <- union(fcres$gene_id[fcres$r1], sares$gene_id[sares$r1])
f <- fcres[match(ids, fcres$gene_id), ]; s <- sares[match(ids, sares$gene_id), ]
t <- data.frame(gene_id = ids, symbol = ifelse(!is.na(f$gene_name), f$gene_name, s$gene_name), gene_biotype = ifelse(!is.na(f$gene_type), f$gene_type, s$gene_type),
  baseMean_fc = f$baseMean, baseMean_salmon = s$baseMean, log2FC_shrunk_fc = f$log2FoldChange_shrunk, log2FC_shrunk_salmon = s$log2FoldChange_shrunk,
  log2FC_unshrunk_fc = f$log2FoldChange, log2FC_unshrunk_salmon = s$log2FoldChange, padj_fc = f$padj, padj_salmon = s$padj, stringsAsFactors = FALSE)
t$rule1_fc <- ifelse(is.na(f$r1), FALSE, f$r1); t$rule1_salmon <- ifelse(is.na(s$r1), FALSE, s$r1)
t$direction_fc <- ifelse(is.na(t$log2FC_shrunk_fc), NA, ifelse(t$log2FC_shrunk_fc > 0, "up", "down"))
t$direction_salmon <- ifelse(is.na(t$log2FC_shrunk_salmon), NA, ifelse(t$log2FC_shrunk_salmon > 0, "up", "down"))
t$rule1_pass <- t$rule1_fc | t$rule1_salmon                                                     # TRUE for every row by construction
t$rule2_pass <- !is.na(t$baseMean_fc) & !is.na(t$baseMean_salmon) & t$baseMean_fc >= BM_MIN & t$baseMean_salmon >= BM_MIN
t$baseMean_ge_50_both <- !is.na(t$baseMean_fc) & !is.na(t$baseMean_salmon) & t$baseMean_fc >= BM_STRICT & t$baseMean_salmon >= BM_STRICT   # stricter tier, not a filter
t$rule3a_pass <- t$rule1_fc & t$rule1_salmon & !is.na(t$direction_fc) & !is.na(t$direction_salmon) & t$direction_fc == t$direction_salmon
lp <- loo_pass[match(ids, rownames(loo_pass)), , drop = FALSE]; tested <- !is.na(match(ids, rownames(loo_pass)))
t$loo_n_pass <- ifelse(tested, rowSums(lp), NA); t$loo_fraction <- t$loo_n_pass / 19; t$loo_n_fail <- 19 - t$loo_n_pass
t$loo_failed_samples <- vapply(seq_along(ids), function(i) if (!tested[i]) "not tested in 19-sample fit" else paste(colnames(lp)[!lp[i, ]], collapse = ";"), "")
t$rule3b_pass <- !is.na(t$loo_n_pass) & t$loo_n_pass >= LOO_MIN_PASS
t$all_rules_pass <- t$rule1_pass & t$rule2_pass & t$rule3a_pass & t$rule3b_pass
why <- function(i) { r <- character(0)
  if (!t$rule2_pass[i]) r <- c(r, sprintf("rule2: baseMean fc %s / salmon %s (need >= %g on both)", format(t$baseMean_fc[i], digits = 3), format(t$baseMean_salmon[i], digits = 3), BM_MIN))
  if (!t$rule3a_pass[i]) r <- c(r, if (!t$rule1_fc[i]) "rule3a: fails rule 1 on featureCounts" else if (!t$rule1_salmon[i]) "rule3a: fails rule 1 on Salmon" else "rule3a: opposite direction")
  if (!t$rule3b_pass[i]) r <- c(r, sprintf("rule3b: passes %s of %d leave-one-out refits (need >= %d)", t$loo_n_pass[i], N_LOO, LOO_MIN_PASS))
  paste(r, collapse = "; ") }
t$fail_reason <- vapply(seq_len(nrow(t)), why, "")
t <- t[order(!t$all_rules_pass, -t$loo_n_pass, t$padj_fc), ]
write.table(t, file.path(out, "candidate_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- funnel (sequential), all genes and protein-coding
pc <- t$gene_biotype == "protein_coding"
step <- function(n, k) data.frame(step = n, all_genes = sum(k), protein_coding = sum(k & pc))
fun <- rbind(step("rule 1 on featureCounts", t$rule1_fc), step("rule 1 on Salmon", t$rule1_salmon), step("rule 1 on either route (candidate table)", t$rule1_pass),
  step(sprintf("+ rule 2 (baseMean >= %g on both routes)", BM_MIN), t$rule1_pass & t$rule2_pass), step("+ rule 3a (rule 1 on both routes, same direction)", t$rule1_pass & t$rule2_pass & t$rule3a_pass),
  step(sprintf("+ rule 3b (LOO >= %d of %d, featureCounts)", LOO_MIN_PASS, N_LOO), t$all_rules_pass), step(sprintf("  of which baseMean >= %g on both (stricter tier)", BM_STRICT), t$all_rules_pass & t$baseMean_ge_50_both))
write.table(fun, file.path(out, "funnel_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(fun)

# ---- per-dropped-sample failures: among featureCounts rule-1 genes (19-sample run), and among the genes that pass rules 1, 2, 3a
fc1 <- fcres$gene_id[fcres$r1]; k <- t$gene_id[t$rule1_pass & t$rule2_pass & t$rule3a_pass]
bys <- data.frame(srr = colnames(cn), label = ss$label, study = as.character(ss$study), group = as.character(ss$group),
  fail_among_fc_rule1_genes = sapply(colnames(cn), function(x) sum(!loo_pass[fc1, x])), n_fc_rule1_genes = length(fc1),
  fail_among_rule123a_genes = sapply(colnames(cn), function(x) sum(!loo_pass[k, x])), n_rule123a_genes = length(k), row.names = NULL)
write.table(bys, file.path(out, "loo_failures_by_dropped_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(bys)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(c(paste("R", getRversion()), paste("DESeq2", packageVersion("DESeq2")), paste("apeglm", packageVersion("apeglm")), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions.txt"))
