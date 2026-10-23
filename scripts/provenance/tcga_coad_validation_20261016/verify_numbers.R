#!/usr/bin/env Rscript
# Independent second implementation (R) of D2, D3 and the hold-out metrics; does not read compare_hosseini.py or score_holdout.R output.
# usage: verify_numbers.R <run1 dir> <prep dir> <09_verification/input dir>
suppressPackageStartupMessages(library(pROC))
a <- commandArgs(trailingOnly = TRUE); run <- a[1]; prep <- a[2]; hin <- a[3]
d <- read.delim(file.path(run, "out/01_de/deseq2_results_all.tsv"), stringsAsFactors = FALSE)
d <- d[!is.na(d$pvalue) & !is.na(d$gene_name) & d$gene_name != "", ]
d <- d[order(-d$baseMean), ]; d <- d[!duplicated(d$gene_name), ]
d$sig <- !is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange_shrunk) >= 1; d$dir <- ifelse(d$log2FoldChange_shrunk > 0, "up", "down")
s3 <- read.delim(file.path(hin, "hosseini_2023_TableS3_genes.tsv"), stringsAsFactors = FALSE); s3$gene <- trimws(s3$gene); s3 <- s3[!duplicated(s3$gene), ]
m <- merge(d, s3, by.x = "gene_name", by.y = "gene"); ms <- m[m$sig, ]
cat(sprintf("D2: significant in S3 %d, same direction %d, agreement %.4f\n", nrow(ms), sum(ms$dir == ms$direction_in_cancer), mean(ms$dir == ms$direction_in_cancer)))
s7 <- strsplit(readLines(file.path(hin, "MOESM2_sheets/Table_S7_p_value_.tsv"))[2], "\t")[[1]]; s7 <- trimws(s7[nzchar(trimws(s7))])
z <- d[match(s7, d$gene_name), ]; test <- !is.na(z$padj); pub <- s3$direction_in_cancer[match(s7, s3$gene)]
ok <- test & z$sig %in% TRUE & z$dir == pub
cat(sprintf("D3: S7 genes %d, testable %d, significant in published direction %d, fraction %.4f; untestable: %s\n", length(s7), sum(test), sum(ok), sum(ok) / sum(test), paste(s7[!test], collapse = ",")))
cat("opposite-direction significant S3 genes:", sum(ms$dir != ms$direction_in_cancer), "\n")
cl <- read.delim(file.path(prep, "clinical_tcga_coad.tsv"), stringsAsFactors = FALSE); sp <- read.delim(file.path(prep, "split.tsv"), stringsAsFactors = FALSE)
cat("patients in both splits:", length(intersect(sp$patient_id[sp$split == "training"], sp$patient_id[sp$split == "holdout"])), "\n")
sc <- read.delim(file.path(run, "../holdout_scoring/holdout_scores.tsv"), stringsAsFactors = FALSE)
y <- as.integer(sc$truth == "tumour"); r <- roc(y, sc$score, levels = c(0, 1), direction = "<", quiet = TRUE)
cat(sprintf("hold-out: n %d, AUC %.6f, min tumour score %.4f, max normal score %.4f, cutoff-free separation: %s\n", length(y), auc(r), min(sc$score[y == 1]), max(sc$score[y == 0]), min(sc$score[y == 1]) > max(sc$score[y == 0])))
