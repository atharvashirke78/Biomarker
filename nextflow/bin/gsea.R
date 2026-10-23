#!/usr/bin/env Rscript
# NEW. usage: gsea.R <deseq2_results_all.tsv> <minSize> <maxSize> <nPermSimple> <padj> <seed> <plot_top_n> [run_kegg true|false, default false]
# run_kegg (params.run_kegg, default false): KEGG licence for commercial use is unresolved. With true the collections and their order are exactly
# as before (Hallmark, KEGG, GO_BP; one set.seed before the loop). With false KEGG is dropped; Hallmark is unaffected, but GO_BP then runs on a
# different position of the random stream, so its permutation-based values may differ slightly from a run with KEGG (see run_log).
suppressPackageStartupMessages({ library(fgsea); library(ggplot2); library(BiocParallel) })
source(file.path(Sys.getenv("PROJECT_BIN", unset = dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), "msig_helper.R"))
a <- commandArgs(trailingOnly = TRUE); mn <- as.integer(a[2]); mx <- as.integer(a[3]); np <- as.integer(a[4]); pc <- as.numeric(a[5]); seed <- as.integer(a[6]); top_n <- as.integer(a[7])
d <- read.delim(a[1], stringsAsFactors = FALSE); d <- d[!is.na(d$stat) & !is.na(d$gene_name) & d$gene_name != "", ]
d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$gene_name), ]                 # one value per symbol: largest |Wald stat|
ranks <- setNames(d$stat, d$gene_name); ranks <- sort(ranks, decreasing = TRUE)
message("ranked genes: ", length(ranks))
run_kegg <- length(a) >= 8 && tolower(a[8]) == "true"
collections <- list(Hallmark = function() get_sets("H"), KEGG = kegg_sets, GO_BP = function() get_sets("C5", "GO:BP"))
if (!run_kegg) { collections$KEGG <- NULL; message("KEGG skipped (--run_kegg false, the default; KEGG licence for commercial use is unresolved)") }
set.seed(seed)
for (nm in names(collections)) {
  ps <- sets_to_list(collections[[nm]]())
  fg <- fgseaMultilevel(ps, ranks, minSize = mn, maxSize = mx, nPermSimple = np, BPPARAM = SerialParam())
  fg <- fg[order(fg$padj), ]; fg$leadingEdge <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
  write.table(as.data.frame(fg), sprintf("gsea_%s.tsv", nm), sep = "\t", quote = FALSE, row.names = FALSE)
  s <- as.data.frame(fg); s <- s[!is.na(s$padj) & s$padj < pc, ]
  message(nm, ": ", length(ps), " sets, ", nrow(s), " with padj < ", pc)
  if (nrow(s) > 0) { s <- s[order(-abs(s$NES)), ]; s <- head(s, top_n); s$pathway <- factor(s$pathway, levels = rev(s$pathway))
    ggsave(sprintf("gsea_%s_top.png", nm), ggplot(s, aes(NES, pathway, fill = NES > 0)) + geom_col() + scale_fill_manual(values = c("blue3", "red3"), guide = "none") +
           theme_bw() + labs(title = paste("fgsea", nm, "(top", top_n, "by |NES|, padj <", pc, ")")), width = 10, height = 7, dpi = 150) }
}
