#!/usr/bin/env Rscript
# usage: cohort_qc.R <counts.tsv|salmon_tximport.rds> <sample_map.tsv (srr,label,study,group)> <sample_qc.tsv> <min_count> <min_samples> <outlier_mad>; sample IDs = SRR
# (refactor Step 3: the three numbers are parameters cohort_min_count / cohort_min_samples / cohort_outlier_mad from conf/thresholds.config)
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(pheatmap); library(jsonlite) })
a <- commandArgs(trailingOnly = TRUE)
min_count <- as.numeric(a[4]); min_samples <- as.numeric(a[5]); mad_cut <- as.numeric(a[6])
ss <- read.delim(a[2], stringsAsFactors = FALSE); rownames(ss) <- ss$srr
if (grepl("\\.rds$", a[1])) { txi <- readRDS(a[1]); cn <- txi$counts
  ss <- ss[colnames(cn), , drop = FALSE]; ss$study <- factor(ss$study); ss$group <- factor(ss$group)
  dds <- DESeqDataSetFromTximport(txi, colData = ss, design = ~ 1)
} else { cn <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
  ss <- ss[colnames(cn), , drop = FALSE]; ss$study <- factor(ss$study); ss$group <- factor(ss$group)
  dds <- DESeqDataSetFromMatrix(round(cn), colData = ss, design = ~ 1) }
dds <- dds[rowSums(counts(dds) >= min_count) >= min(min_samples, ncol(dds)), ]
vsd <- tryCatch(vst(dds, blind = TRUE), error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
m <- assay(vsd)
pc <- prcomp(t(m)); pv <- 100 * pc$sdev^2 / sum(pc$sdev^2)
sc <- data.frame(srr = rownames(pc$x), label = ss$label, PC1 = pc$x[, 1], PC2 = pc$x[, 2], study = ss$study, group = ss$group)
ggsave("pca_group_study.png", ggplot(sc, aes(PC1, PC2, colour = group, shape = study, label = srr)) + geom_point(size = 3.5) +
  geom_text(vjust = -1, size = 3, colour = "black") + theme_bw() + xlab(sprintf("PC1: %.1f%%", pv[1])) + ylab(sprintf("PC2: %.1f%%", pv[2])), width = 7, height = 5.5, dpi = 150)
pheatmap(cor(m), annotation_col = ss[, c("study", "group")], filename = "sample_correlation_heatmap.png", width = 8, height = 7)
ls <- data.frame(srr = colnames(cn), lib = colSums(cn), group = ss$group)
ggsave("library_size_barplot.png", ggplot(ls, aes(srr, lib / 1e6, fill = group)) + geom_col() + theme_bw() + ylab("assigned counts (millions)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)), width = 8, height = 4.5, dpi = 150)
# flag only: > 3 MAD from the group centre (median) on PC1 or PC2, computed WITHIN EACH STUDY (study x group cell), so the
# study (batch) effect on PC1 does not flag a whole study. Nothing is removed.
rz <- function(x) { s <- mad(x); if (s == 0) s <- sd(x); if (is.na(s) || s == 0) s <- 1; (x - median(x)) / s }
cell <- interaction(sc$study, sc$group, drop = TRUE)
sc$z_pc1 <- ave(sc$PC1, cell, FUN = rz); sc$z_pc2 <- ave(sc$PC2, cell, FUN = rz)
sc$flag_outlier_3mad <- abs(sc$z_pc1) > mad_cut | abs(sc$z_pc2) > mad_cut   # column name kept (contract qc_table reads it); the limit is cohort_outlier_mad
sc$lib_size <- ls$lib[match(sc$srr, ls$srr)]
sq <- read.delim(a[3], stringsAsFactors = FALSE); sc <- merge(sc, sq[, c("srr", "status")], by = "srr", all.x = TRUE)
# ---- NEW (R7): variance partition on PC1/PC2: R-squared of a one-way ANOVA per PC for study (batch) and for group. Report/flag only.
r2_one <- function(y, f) { f <- droplevels(factor(f)); if (nlevels(f) < 2) return(NA_real_)
  1 - sum(tapply(y, f, function(v) sum((v - mean(v))^2))) / sum((y - mean(y))^2) }
vp <- lapply(1:min(2, ncol(pc$x)), function(i) list(study_r2 = r2_one(pc$x[, i], ss$study), group_r2 = r2_one(pc$x[, i], ss$group)))
names(vp) <- paste0("PC", seq_along(vp))
n_study_levels <- nlevels(droplevels(ss$study))
study_dominates_group <- isTRUE(n_study_levels >= 2 && !is.na(vp$PC1$study_r2) && (is.na(vp$PC1$group_r2) || vp$PC1$study_r2 > vp$PC1$group_r2))
sc$study_dominates_group <- study_dominates_group     # cohort-level flag repeated on every row (new last column)
write.table(sc, "cohort_qc_flags.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
r2 <- sapply(1:min(4, ncol(pc$x)), function(i) c(study = summary(lm(pc$x[, i] ~ ss$study))$r.squared, group = summary(lm(pc$x[, i] ~ ss$group))$r.squared))
res <- list(n_samples = ncol(m), n_genes = nrow(m), percent_var = round(pv[1:min(4, length(pv))], 1), r2_by_pc = as.data.frame(t(r2)),
                outliers_flagged = I(sc$srr[sc$flag_outlier_3mad]), flag_rule = sprintf("|robust z| > %g (MAD) on PC1 or PC2 within each study x group cell", mad_cut), min_count = min_count, min_samples = min_samples, outlier_mad = mad_cut, note = "flag only, nothing removed")

# ---- NEW (H4b): sex check from the count matrix, flag only, never excludes. XIST and Y-chromosome genes by symbol via the pipeline symbol table.
sym_file <- Sys.getenv("COHORT_QC_SYMBOL_TABLE", "/workspace/Biomarker_pipeline/ref/gene_id_to_symbol.tsv")
ychr_genes <- c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D")
sex_note <- character(0)
l2cpm <- function(ids) { ids <- ids[ids %in% rownames(cn)]; if (!length(ids)) return(NULL)
  log2(sweep(cn[ids, , drop = FALSE], 2, colSums(cn), "/") * 1e6 + 1) }
xist_v <- rep(NA_real_, ncol(cn)); y_v <- rep(NA_real_, ncol(cn)); genes_found <- character(0)
if (file.exists(sym_file)) {
  st <- tryCatch(read.delim(sym_file, stringsAsFactors = FALSE, quote = ""), error = function(e) NULL)
  if (!is.null(st) && all(c("gene_id", "gene_name") %in% colnames(st))) {
    ids_of <- function(g) { x <- st$gene_id[st$gene_name %in% g]; unique(c(x, intersect(rownames(cn), sub("\\.[0-9]+$", "", x)))) }
    xm <- l2cpm(ids_of("XIST")); if (!is.null(xm)) { xist_v <- colMeans(xm); genes_found <- c(genes_found, "XIST") }
    for (g in ychr_genes) if (!is.null(l2cpm(ids_of(g)))) genes_found <- c(genes_found, g)
    ym <- l2cpm(ids_of(ychr_genes)); if (!is.null(ym)) y_v <- colMeans(ym)
  } else sex_note <- c(sex_note, "symbol table unreadable or lacks gene_id/gene_name")
} else sex_note <- c(sex_note, paste("symbol table not found:", sym_file))
if (!"XIST" %in% genes_found) sex_note <- c(sex_note, "XIST not in the count matrix")
if (!any(ychr_genes %in% genes_found)) sex_note <- c(sex_note, "no Y-chromosome gene from the list in the count matrix")
# call only if the score (XIST - mean Y log2CPM; XIST alone if no Y gene) is clearly bimodal: the largest gap between sorted values must be
# >= 3 log2 units and more than twice the next-largest gap (a real separation, not one step in a continuum), each cluster >= 1 sample. Otherwise 'not called'.
score <- if (all(is.na(y_v))) xist_v else xist_v - y_v
call <- rep("not called", length(score)); sex_rule <- "not called: fewer than 3 samples or no usable score"
if (sum(!is.na(score)) >= 3) {
  o <- order(score); sv <- score[o]; gaps <- diff(sv); k <- which.max(gaps); lo <- sv[1:k]; hi <- sv[(k + 1):length(sv)]
  second <- if (length(gaps) > 1) max(gaps[-k]) else 0
  if (gaps[k] >= 3 && gaps[k] > 2 * second) { call[score <= sv[k]] <- "male"; call[score >= sv[k + 1]] <- "female"
    sex_rule <- sprintf("called: bimodal score, gap %.2f log2 units between clusters (low XIST/high Y = male)", gaps[k])
  } else sex_rule <- sprintf("not called: score not clearly bimodal (largest gap %.2f, needs >= 3 and > 2 x the next-largest gap %.2f)", gaps[k], second)
}
sexcol <- intersect(c("sex", "Sex", "SEX", "gender"), c(colnames(ss), colnames(sq)))
meta_sex <- rep(NA_character_, ncol(cn))
if (length(sexcol)) { v <- if (sexcol[1] %in% colnames(ss)) ss[colnames(cn), sexcol[1]] else sq[match(colnames(cn), sq$srr), sexcol[1]]
  v <- tolower(trimws(as.character(v))); meta_sex <- ifelse(v %in% c("f", "female", "woman"), "female", ifelse(v %in% c("m", "male", "man"), "male", NA)) }
mismatch <- ifelse(is.na(meta_sex) | call == "not called", NA, meta_sex != call)
sexdf <- data.frame(srr = colnames(cn), label = ss$label, xist_log2cpm = round(xist_v, 3), ychr_mean_log2cpm = round(y_v, 3), sex_call = call,
                    sex_metadata = meta_sex, sex_mismatch_flag = mismatch, stringsAsFactors = FALSE)
write.table(sexdf, "cohort_qc_sex.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
res$variance_partition_pc <- vp
res$study_levels <- n_study_levels
res$study_dominates_group <- study_dominates_group
res$variance_partition_note <- "R-squared of a one-way ANOVA of each PC score on study and on group; study_dominates_group = study R2 > group R2 on PC1 and >= 2 study levels; flag only"
res$sex_check <- list(genes_found = genes_found, rule = sex_rule, n_called = sum(call != "not called"), sex_column_used = if (length(sexcol)) sexcol[1] else NA,
                      mismatch_samples = I(sexdf$srr[which(sexdf$sex_mismatch_flag)]), notes = I(sex_note), note = "flag only, nothing removed")
write_json(res, "cohort_qc_summary.json", auto_unbox = TRUE, pretty = TRUE, na = "null")
