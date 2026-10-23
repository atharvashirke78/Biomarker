#!/usr/bin/env Rscript
# Step 1: per-sample QC table. Values come from the primary files MultiQC summarises (STAR Log.final.out, featureCounts .summary,
# Salmon meta_info.json / lib_format_counts.json, fastp json) and are cross-checked against results/04_counts/sample_qc.tsv.
# Cutoffs are read from nextflow/nextflow.config (qc_*), never changed.
suppressPackageStartupMessages(library(jsonlite))
B <- "/workspace/Biomarker_pipeline"; R <- file.path(B, "results"); O <- file.path(R, "09_verification/step1_qc")
sm <- read.delim(file.path(R, "sample_map.tsv"), stringsAsFactors = FALSE)
sq <- read.delim(file.path(R, "04_counts/sample_qc.tsv"), stringsAsFactors = FALSE); rownames(sq) <- sq$srr
cfg <- readLines(file.path(B, "nextflow/nextflow.config"))
cut <- function(k) as.numeric(sub("^.*=\\s*([0-9.]+).*$", "\\1", grep(paste0("^\\s*", k, "\\s*="), cfg, value = TRUE)[1]))
cs <- c(star_warn = cut("qc_star_warn"), star_fail = cut("qc_star_fail"), fc_warn = cut("qc_fc_warn"), salmon_warn = cut("qc_salmon_warn"), salmon_fail = cut("qc_salmon_fail"))
rows <- lapply(sm$srr, function(s) {
  lg <- readLines(file.path(R, "02_alignment", s, "Log.final.out"))
  star <- as.numeric(sub("%", "", sub(".*\\|\\s*", "", grep("Uniquely mapped reads %", lg, value = TRUE))))
  inreads <- as.numeric(sub(".*\\|\\s*", "", grep("Number of input reads", lg, value = TRUE)))
  fs <- read.delim(list.files(file.path(R, "04_counts", s), pattern = "summary$", full.names = TRUE), stringsAsFactors = FALSE)
  fc_assigned <- fs[fs[[1]] == "Assigned", 2]; fc <- 100 * fc_assigned / sum(fs[[2]])
  mi <- fromJSON(file.path(R, "03_quant_salmon", s, "aux_info/meta_info.json"))
  lf <- fromJSON(file.path(R, "03_quant_salmon", s, "lib_format_counts.json"))
  fp <- fromJSON(file.path(R, "01_fastp", s, "fastp.json"))
  isr <- lf$ISR; isf <- lf$ISF; ratio <- if ((isr + isf) > 0) isr / (isr + isf) else NA
  strand <- if (is.na(ratio)) "NA" else if (ratio >= 0.9) "reverse (ISR)" else if (ratio <= 0.1) "forward (ISF)" else "unstranded"
  data.frame(srr = s, label = sm$label[sm$srr == s], study = sm$study[sm$srr == s], group = sm$group[sm$srr == s],
             input_pairs_M = round(inreads / 1e6, 2), star_unique_pct = star, fc_assigned_pct = round(fc, 2), salmon_mapped_pct = round(mi$percent_mapped, 2),
             fastp_q30_pct = round(100 * fp$summary$after_filtering$q30_rate, 2), library_size_M = round(fc_assigned / 1e6, 2),
             fastp_duplication_pct = round(100 * fp$duplication$rate, 1), salmon_ISR_frac = round(ratio, 3), inferred_strandedness = strand,
             salmon_expected_format = lf$expected_format, stringsAsFactors = FALSE)
})
t <- do.call(rbind, rows)
# consistency with the pipeline's own sample_qc.tsv
t$matches_sample_qc <- abs(t$star_unique_pct - sq[t$srr, "star_unique_pct"]) < 0.01 & abs(t$fc_assigned_pct - sq[t$srr, "fc_assigned_pct"]) < 0.01 &
                       abs(t$salmon_mapped_pct - sq[t$srr, "salmon_mapped_pct"]) < 0.01 & abs(t$fastp_q30_pct - sq[t$srr, "fastp_q30_pct"]) < 0.01
t$pipeline_status <- sq[t$srr, "status"]
why <- function(i) { w <- c()
  if (t$star_unique_pct[i] < cs["star_fail"]) w <- c(w, sprintf("STAR unique %.2f < %g (FAIL)", t$star_unique_pct[i], cs["star_fail"])) else if (t$star_unique_pct[i] < cs["star_warn"]) w <- c(w, sprintf("STAR unique %.2f < %g (WARN)", t$star_unique_pct[i], cs["star_warn"]))
  if (t$fc_assigned_pct[i] < cs["fc_warn"]) w <- c(w, sprintf("featureCounts assigned %.2f < %g (WARN)", t$fc_assigned_pct[i], cs["fc_warn"]))
  if (t$salmon_mapped_pct[i] < cs["salmon_fail"]) w <- c(w, sprintf("Salmon mapped %.2f < %g (FAIL)", t$salmon_mapped_pct[i], cs["salmon_fail"])) else if (t$salmon_mapped_pct[i] < cs["salmon_warn"]) w <- c(w, sprintf("Salmon mapped %.2f < %g (WARN)", t$salmon_mapped_pct[i], cs["salmon_warn"]))
  paste(w, collapse = "; ") }
t$outside_cutoffs <- vapply(seq_len(nrow(t)), why, "")
write.table(t, file.path(O, "per_sample_qc_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
out <- t[t$outside_cutoffs != "", c("srr", "label", "study", "group", "pipeline_status", "outside_cutoffs")]
txt <- c(sprintf("cutoffs (nextflow.config): %s", paste(names(cs), cs, sep = "=", collapse = ", ")),
         sprintf("samples: %d; table values equal the pipeline's sample_qc.tsv for %d of %d", nrow(t), sum(t$matches_sample_qc), nrow(t)),
         sprintf("samples outside a cutoff: %d", nrow(out)), capture.output(print(out, row.names = FALSE, right = FALSE)), "",
         "inferred strandedness (Salmon ISR/(ISR+ISF)): ", capture.output(print(table(t$inferred_strandedness))),
         sprintf("fastp duplication %% per study (median): %s", paste(names(tapply(t$fastp_duplication_pct, t$study, median)), tapply(t$fastp_duplication_pct, t$study, median), sep = "=", collapse = ", ")),
         sprintf("library size M per study (median): %s", paste(names(tapply(t$library_size_M, t$study, median)), tapply(t$library_size_M, t$study, median), sep = "=", collapse = ", ")))
writeLines(txt, file.path(O, "qc_outside_cutoffs.txt")); cat(txt, sep = "\n")
