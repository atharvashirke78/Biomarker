# Compares the new deseq2_design.R output (design ~ study + group, 19 samples) with the main-run results.
a <- read.delim("results/06_dge/deseq2_results_all.tsv", stringsAsFactors = FALSE)
b <- read.delim("results/10_independent_validation/script_check_main19/main19_rerun_results_all.tsv", stringsAsFactors = FALSE)
cat("rows main", nrow(a), "new", nrow(b), "; same gene_id set:", setequal(a$gene_id, b$gene_id), "\n")
m <- merge(a, b, by = "gene_id", suffixes = c(".a", ".b"))
for (k in c("baseMean", "log2FoldChange", "log2FoldChange_shrunk", "pvalue", "padj")) {
  x <- m[[paste0(k, ".a")]]; y <- m[[paste0(k, ".b")]]
  cat(sprintf("%-22s NA-pattern same: %s; max abs diff: %g\n", k, identical(is.na(x), is.na(y)), max(abs(x - y), na.rm = TRUE)))
}
sg <- function(d) d$gene_id[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1]
cat("significant: main", length(sg(a)), "new", length(sg(b)), "identical sets:", setequal(sg(a), sg(b)), "\n")
