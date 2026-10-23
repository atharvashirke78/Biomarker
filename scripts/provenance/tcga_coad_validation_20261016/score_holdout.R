#!/usr/bin/env Rscript
# PREREG (b)/(e)/D4: score the hold-out ONCE with the model genes, signs and cutoff of the training run (no refit, no re-tuning).
# Score = the formula in model_spec.txt (same as fit_diagnostic_model.R sc_fun): mean log2(CPM+1) of UP genes - mean of DOWN genes,
# CPM total = the sample's column sum over the genes of the DISCOVER input matrix (the same rows as counts_training.tsv).
# usage: score_holdout.R <counts.tsv> <clinical.tsv> <holdout_samples.txt> <run 03_model dir> <out dir (must not exist)>
suppressPackageStartupMessages(library(pROC))
a <- commandArgs(trailingOnly = TRUE); out <- a[5]; if (dir.exists(out)) stop(out, " exists: refusing to overwrite"); dir.create(out)
cn <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
cl <- read.delim(a[2], stringsAsFactors = FALSE, colClasses = "character"); rownames(cl) <- cl$sample_id
ho <- readLines(a[3]); ho <- ho[nzchar(ho)]; stopifnot(all(ho %in% colnames(cn)), all(ho %in% cl$sample_id))
mg <- read.delim(file.path(a[4], "model_genes.tsv"), stringsAsFactors = FALSE)
spec <- readLines(file.path(a[4], "model_spec.txt"))
cut_line <- grep("^cutoff = ", spec, value = TRUE); stopifnot(length(cut_line) == 1)
cutoff <- as.numeric(sub("^cutoff = ([^ ]+).*$", "\\1", cut_line)); stopifnot(is.finite(cutoff))
stopifnot(all(mg$gene_id %in% rownames(cn)))
up <- mg$gene_id[mg$set == "up"]; dn <- mg$gene_id[mg$set == "down"]
L4 <- log2(sweep(cn[, ho, drop = FALSE], 2, colSums(cn[, ho, drop = FALSE]), "/") * 1e6 + 1)
score <- (if (length(up)) colMeans(L4[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(L4[dn, , drop = FALSE]) else 0)
y <- as.integer(cl[ho, "group"] == "tumour"); pred <- as.integer(score > cutoff)
r <- roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE); ci <- ci.auc(r, method = "delong")
TP <- sum(pred == 1 & y == 1); FN <- sum(pred == 0 & y == 1); TN <- sum(pred == 0 & y == 0); FP <- sum(pred == 1 & y == 0)
se <- binom.test(TP, TP + FN)$conf.int; sp <- binom.test(TN, TN + FP)$conf.int
P <- data.frame(split = "HOLD-OUT (scored once; model and cutoff from training)", n = length(y), n_tumour = sum(y), n_normal = sum(1 - y),
  AUC = as.numeric(auc(r)), AUC_delong_ci_low = ci[1], AUC_delong_ci_high = ci[3], cutoff = cutoff, TP = TP, FN = FN, TN = TN, FP = FP,
  sensitivity = TP / (TP + FN), sens_cp_ci_low = se[1], sens_cp_ci_high = se[2], specificity = TN / (TN + FP), spec_cp_ci_low = sp[1], spec_cp_ci_high = sp[2],
  balanced_accuracy = (TP / (TP + FN) + TN / (TN + FP)) / 2, n_up_genes = length(up), n_down_genes = length(dn))
write.table(P, file.path(out, "holdout_performance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample_id = ho, patient_id = cl[ho, "patient_id"], truth = cl[ho, "group"], score = score, pred_tumour = pred),
            file.path(out, "holdout_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("R", getRversion()), paste("pROC", packageVersion("pROC"))), file.path(out, "tool_versions.txt"))
print(t(P))
