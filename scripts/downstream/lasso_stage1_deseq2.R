# rc1-prep copy (tracked) of results/14_lasso_trial/lasso_stage1_deseq2.R; original left unchanged as the record of the old run. Args: no defaults; refuses a non-empty output folder unless --force.
#!/usr/bin/env Rscript
# LASSO trial, stage 1 (needs DESeq2: env/prefix/r_dge). Run from the project root. featureCounts route, 19 samples, seed 42.
# For every fold (19 leave-one-out + 2 leave-one-study-out) it (1) refits DESeq2 on the TRAINING samples only (same call as results/11_candidate_selection:
# ~ study + group, or ~ group for a one-study training set; pre-filter >= 10 counts in >= n samples (n = 9, full cohort, fixed); results(alpha = 0.05); lfcShrink apeglm),
# (2) keeps genes with rule 1 (padj < 0.05, |shrunk log2FC| >= 1) and rule 2 (baseMean >= 20) from THAT refit, (3) makes the two normalisation arms.
# Arm 1 (leaky, comparison only): one VST (blind = TRUE) over all 19 samples.
# Arm 2: VST fitted on the training samples only (blind = TRUE). Test samples: size factor = median of ratios against the TRAINING geometric means
#   (computed directly as DESeq2 does: exp(median(log(count) - log(training geomean)) over genes with finite geomean and count > 0; DESeq2's own estimateSizeFactors(geoMeans=) returns 1 for a single-sample
#   object, which was caught by the sanity check, so the formula is applied by hand and set with sizeFactors<-), then the training dispersion function is frozen (dispersionFunction(dds_test) <- dispersionFunction(dds_train))
#   and varianceStabilizingTransformation(blind = FALSE) applies the training parametric VST formula (no re-estimation on test data).
# Run: env/prefix/r_dge/bin/Rscript scripts/downstream/lasso_stage1_deseq2.R <counts_dir> <out_dir> [--force]. Reads only <counts_dir>/counts_featurecounts.tsv and sample_map.tsv (parent of counts_dir); alpha, count and n from conf/thresholds.config. The 15 validation samples are not used.
suppressPackageStartupMessages({ library(DESeq2); library(apeglm) })
set.seed(42)
args <- commandArgs(trailingOnly = TRUE); force <- "--force" %in% args; args <- setdiff(args, "--force")
usage <- "usage: %s %s [--force]   (no defaults)"
if (length(args) != 2) stop(sprintf(usage, "lasso_stage1_deseq2.R", "<counts_dir> <out_dir>"))
counts_dir <- args[1]; out <- args[2]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) > 0 && !force) stop("ERROR: output folder ", out, " already contains files; nothing written. Use a new folder (or --force).")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (d in c("refits", "features", "arm2")) dir.create(file.path(out, d), showWarnings = FALSE)
yaml_vals <- function(file, block, key) { l <- readLines(file); i <- grep(paste0("^", block, ":"), l); stopifnot(length(i) == 1)
  j <- i + 1; while (j <= length(l) && grepl("^(\\s|#|$)", l[j])) j <- j + 1; blk <- l[(i + 1):(j - 1)]
  v <- grep(paste0("^\\s*", key, ":"), blk, value = TRUE); as.numeric(sub("^\\s*[a-z0-9_]+:\\s*([0-9.]+).*$", "\\1", v)) }   # numbers of a key inside a top-level yaml block (in file order)
RULES <- "config/candidate_rules.yaml"                                           # rules 1-2 are read from here, not typed
P_CUT <- yaml_vals(RULES, "rule1_significance", "padj_max"); LFC_CUT <- yaml_vals(RULES, "rule1_significance", "abs_shrunk_log2fc_min")
BM_MIN <- yaml_vals(RULES, "rule2_measurable", "basemean_min")[1]
stopifnot(length(P_CUT) == 1, length(LFC_CUT) == 1, is.finite(c(P_CUT, LFC_CUT, BM_MIN)))
read_key <- function(file, key, sep) { l <- grep(paste0("^\\s*", key, "\\s*", sep), readLines(file), value = TRUE); stopifnot(length(l) == 1)
  trimws(gsub("['\"]", "", sub("(#|//).*$", "", sub(paste0("^\\s*", key, "\\s*", sep), "", l)))) }   # same helper as scripts/candidate_selection.R
ALPHA <- as.numeric(read_key("conf/thresholds.config", "de_padj", "=")); MIN_COUNT <- as.numeric(read_key("conf/thresholds.config", "de_min_count", "=")); DE_MIN_SAMPLES <- read_key("conf/thresholds.config", "de_min_samples", "=")
stopifnot(is.finite(ALPHA), is.finite(MIN_COUNT))
cn <- as.matrix(read.delim(file.path(counts_dir, "counts_featurecounts.tsv"), row.names = 1, check.names = FALSE)); stopifnot(ncol(cn) == 19)
sm <- read.delim(file.path(dirname(normalizePath(counts_dir)), "sample_map.tsv"), stringsAsFactors = FALSE); rownames(sm) <- sm$srr
ss <- sm[colnames(cn), ]; ss$study <- factor(ss$study); ss$group <- relevel(factor(ss$group), ref = "normal")
MIN_SAMPLES <- if (DE_MIN_SAMPLES == "auto") min(table(ss$group)) else as.numeric(DE_MIN_SAMPLES)   # from the FULL 19-sample cohort, fixed for every fold (not recalculated per training set)
cat(sprintf("settings: alpha %g | pre-filter >= %g counts in >= %g samples (fixed for all folds)\n", ALPHA, MIN_COUNT, MIN_SAMPLES))
cn <- round(cn)
mk <- function(cols, design) { d <- DESeqDataSetFromMatrix(cn[, cols, drop = FALSE], colData = ss[cols, ], design = design); d$study <- droplevels(d$study); d$group <- droplevels(d$group); d }
prefilter <- function(d) d[rowSums(counts(d) >= MIN_COUNT) >= MIN_SAMPLES, ]

sf_ref <- function(m, gm) { lg <- log(gm); apply(m, 2, function(x) exp(median((log(x) - lg)[is.finite(lg) & x > 0]))) }

# ---- Arm 1: one VST over all 19 samples (leaky)
d19 <- prefilter(mk(colnames(cn), ~ study + group)); d19 <- estimateSizeFactors(d19)
v19 <- varianceStabilizingTransformation(d19, blind = TRUE)
write.table(data.frame(gene_id = rownames(d19), assay(v19), check.names = FALSE), file.path(out, "arm1_vst_all19.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

folds <- c(setNames(lapply(colnames(cn), function(h) list(train = setdiff(colnames(cn), h), test = h)), paste0("loo_", colnames(cn))),
  los_A_to_B = list(list(train = rownames(ss)[ss$study == "SRP219837"], test = rownames(ss)[ss$study == "SRP301216"])),
  los_B_to_A = list(list(train = rownames(ss)[ss$study == "SRP301216"], test = rownames(ss)[ss$study == "SRP219837"])))
man <- list()
for (fn in names(folds)) {
  tr <- folds[[fn]]$train; te <- folds[[fn]]$test; message(fn, ": train ", length(tr), ", test ", length(te))
  design <- if (nlevels(droplevels(ss[tr, "study"])) > 1) ~ study + group else ~ group
  d <- prefilter(mk(tr, design)); d <- DESeq(d, quiet = TRUE)
  res <- results(d, name = "group_tumour_vs_normal", alpha = ALPHA); shr <- lfcShrink(d, coef = "group_tumour_vs_normal", type = "apeglm", res = res, quiet = TRUE)
  gt <- data.frame(gene_id = rownames(res), baseMean = res$baseMean, lfc = res$log2FoldChange, lfc_shr = shr$log2FoldChange, padj = res$padj, stringsAsFactors = FALSE)
  gz <- gzfile(file.path(out, "refits", paste0(fn, "_gene_table.tsv.gz")), "w"); write.table(gt, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  feat <- gt$gene_id[!is.na(gt$padj) & gt$padj < P_CUT & !is.na(gt$lfc_shr) & abs(gt$lfc_shr) >= LFC_CUT & gt$baseMean >= BM_MIN]
  writeLines(feat, file.path(out, "features", paste0(fn, ".txt")))
  # Arm 2: training-only VST (design ~1, blind) on the same pre-filtered training genes
  dt <- prefilter(mk(tr, ~ 1)); dt <- estimateSizeFactors(dt); dt <- estimateDispersionsGeneEst(dt, quiet = TRUE); dt <- estimateDispersionsFit(dt, quiet = TRUE, fitType = "parametric")
  ft <- attr(dispersionFunction(dt), "fitType")
  if (ft != "parametric") {   # DESeq2 fell back to a local fit (too few training samples); the local VST cannot be frozen exactly, so use the exact-frozen "mean" fit for this fold and record it
    message(fn, ": parametric dispersion fit not available (", ft, "); using fitType = 'mean' for the frozen arm-2 transform")
    dt <- estimateDispersionsFit(dt, quiet = TRUE, fitType = "mean"); ft <- attr(dispersionFunction(dt), "fitType"); stopifnot(ft == "mean") }
  vtr <- getVarianceStabilizedData(dt)
  gm <- exp(rowMeans(log(counts(dt))))                                         # training geometric means (genes with a zero count get 0 and are ignored by DESeq2)
  dte <- DESeqDataSetFromMatrix(cn[rownames(dt), te, drop = FALSE], colData = ss[te, , drop = FALSE], design = ~ 1)
  sizeFactors(dte) <- sf_ref(counts(dte), gm); dispersionFunction(dte) <- dispersionFunction(dt)
  vte <- assay(varianceStabilizingTransformation(dte, blind = FALSE))
  s1 <- tr[1]; dchk <- DESeqDataSetFromMatrix(cn[rownames(dt), s1, drop = FALSE], colData = ss[s1, , drop = FALSE], design = ~ 1)   # sanity (every fold): a TRAINING sample through the frozen transform must reproduce the training VST
  sizeFactors(dchk) <- sf_ref(counts(dchk), gm); dispersionFunction(dchk) <- dispersionFunction(dt)
  dif <- max(abs(assay(varianceStabilizingTransformation(dchk, blind = FALSE))[, 1] - vtr[, s1])); message("sanity ", fn, " (training sample through frozen transform): max abs diff ", signif(dif, 3)); stopifnot(dif < 1e-8)
  f2 <- intersect(feat, rownames(vtr))
  write.table(data.frame(gene_id = f2, vtr[f2, , drop = FALSE], check.names = FALSE), file.path(out, "arm2", paste0(fn, "_train.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(gene_id = f2, vte[f2, , drop = FALSE], check.names = FALSE), file.path(out, "arm2", paste0(fn, "_test.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  man[[fn]] <- data.frame(fold = fn, n_train = length(tr), n_test = length(te), train = paste(tr, collapse = ";"), test = paste(te, collapse = ";"), design = deparse(design),
    n_genes_after_prefilter = nrow(d), n_rule1_rule2_features = length(feat), n_features_in_arm2_matrix = length(f2), disp_fit_type = ft,
    test_size_factors = paste(signif(sizeFactors(dte), 4), collapse = ";"), train_size_factor_range = paste(signif(range(sizeFactors(dt)), 4), collapse = "-"))
}
write.table(do.call(rbind, man), file.path(out, "fold_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("R", getRversion()), paste("DESeq2", packageVersion("DESeq2")), paste("apeglm", packageVersion("apeglm")), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions_stage1.txt"))
message("stage 1 done")
