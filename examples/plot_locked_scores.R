#!/usr/bin/env Rscript
# Bar chart of the LOCKED-model scores of one apply run against the locked cutoff (read from lock/diagnostic_v1/final_model_spec.txt).
# usage: examples/plot_locked_scores.R <scores.tsv from 09_locked_model/scores> <out.png> [title]
suppressPackageStartupMessages(library(ggplot2))
a <- commandArgs(trailingOnly = TRUE)
sc <- read.delim(a[1], stringsAsFactors = FALSE)
spec <- readLines(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "..", "lock", "diagnostic_v1", "final_model_spec.txt"))
cl <- sub("^cutoff[^0-9-]*", "", grep("^cutoff", spec, value = TRUE)[1]); cut <- as.numeric(regmatches(cl, regexpr("-?[0-9]+\\.?[0-9]*", cl)))
if (length(cut) != 1 || is.na(cut)) stop("cutoff not found in ", "lock/diagnostic_v1/final_model_spec.txt")
sc$sample <- factor(sc$sample, levels = sc$sample[order(sc$score)])
sc$call <- ifelse(sc$call_tumour == 1, "called tumour", "called normal")
p <- ggplot(sc, aes(sample, score, fill = truth)) + geom_col(width = .7) + geom_hline(yintercept = cut, linetype = 2) +
  annotate("text", x = 1, y = cut, label = sprintf("locked cutoff %.3f", cut), vjust = -.6, hjust = 0, size = 3.2) +
  scale_fill_manual(values = c(normal = "#4C78A8", tumour = "#E45756"), name = "true group") + coord_flip() + theme_bw() +
  labs(x = NULL, y = "locked-model score (mean log2 CPM+1 of up genes minus down genes)", title = if (length(a) > 2) a[3] else "Locked-model scores") + theme(plot.title = element_text(size = 10))
ggsave(a[2], p, width = 7.5, height = 3.6, dpi = 150)
