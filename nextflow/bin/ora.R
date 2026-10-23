#!/usr/bin/env Rscript
# NEW. usage: ora.R <deseq2_results_all.tsv> <padj> <lfc> <minGS> <maxGS> <qcut> <dotplot_top_n> [run_kegg true|false, default false]
# run_kegg (params.run_kegg, default false): KEGG licence for commercial use is unresolved, so KEGG ORA is skipped unless explicitly switched on.
# With run_kegg true the code path is exactly the previous one; with false only the KEGG files are not written (ora_kegg_source.txt records the skip).
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2) })
source(file.path(Sys.getenv("PROJECT_BIN", unset = dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), "msig_helper.R"))
a <- commandArgs(trailingOnly = TRUE); pc <- as.numeric(a[2]); lc <- as.numeric(a[3]); mn <- as.integer(a[4]); mx <- as.integer(a[5]); qc <- as.numeric(a[6]); top_n <- as.integer(a[7])
run_kegg <- length(a) >= 8 && tolower(a[8]) == "true"
d <- read.delim(a[1], stringsAsFactors = FALSE)
d <- d[!is.na(d$gene_name) & d$gene_name != "" & !is.na(d$pvalue), ]            # universe: genes tested by DESeq2 after pre-filter
sig <- !is.na(d$padj) & d$padj < pc & abs(d$log2FoldChange) >= lc
sets <- list(up = unique(d$gene_name[sig & d$log2FoldChange > 0]), down = unique(d$gene_name[sig & d$log2FoldChange < 0]))
universe <- unique(d$gene_name)
map <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
uni_e <- unique(map$ENTREZID)
empty <- function(f) write.table(data.frame(ID = character(), Description = character(), p.adjust = numeric()), f, sep = "\t", quote = FALSE, row.names = FALSE)
kegg_src <- "KEGG REST web service (clusterProfiler::enrichKEGG, hsa)"
for (dir in names(sets)) {
  g <- sets[[dir]]; ge <- unique(map$ENTREZID[map$SYMBOL %in% g]); message(dir, ": ", length(g), " symbols, ", length(ge), " Entrez ids")
  fgo <- sprintf("ora_%s_GO_BP.tsv", dir); fk <- sprintf("ora_%s_KEGG.tsv", dir)
  if (length(ge) < 3) { empty(fgo); if (run_kegg) empty(fk); next }
  go <- enrichGO(ge, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = uni_e, pAdjustMethod = "BH", pvalueCutoff = qc, qvalueCutoff = qc, minGSSize = mn, maxGSSize = mx, readable = TRUE)
  if (!is.null(go) && nrow(as.data.frame(go)) > 0) { write.table(as.data.frame(go), fgo, sep = "\t", quote = FALSE, row.names = FALSE)
    ggsave(sprintf("ora_%s_GO_BP_dotplot.png", dir), dotplot(go, showCategory = top_n) + ggtitle(paste("GO BP ORA:", dir)), width = 9, height = 7, dpi = 150) } else empty(fgo)
  if (!run_kegg) next
  kk <- tryCatch(enrichKEGG(ge, organism = "hsa", universe = uni_e, pAdjustMethod = "BH", pvalueCutoff = qc, qvalueCutoff = qc, minGSSize = mn, maxGSSize = mx), error = function(e) { message("enrichKEGG failed: ", conditionMessage(e)); NULL })
  if (is.null(kk)) {   # KEGG web service unreachable -> documented fallback, not a workaround: MSigDB KEGG sets via msigdbr + enricher
    kegg_src <- "FALLBACK: MSigDB C2 CP:KEGG sets via msigdbr + clusterProfiler::enricher (KEGG web service not reachable)"
    ks <- kegg_sets(); t2g <- data.frame(term = ks$gs_name, gene = ks$gene_symbol)
    kk <- enricher(g, universe = universe, TERM2GENE = t2g, pAdjustMethod = "BH", pvalueCutoff = qc, qvalueCutoff = qc, minGSSize = mn, maxGSSize = mx)
  } else kk <- setReadable(kk, org.Hs.eg.db, keyType = "ENTREZID")
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) { write.table(as.data.frame(kk), fk, sep = "\t", quote = FALSE, row.names = FALSE)
    ggsave(sprintf("ora_%s_KEGG_dotplot.png", dir), dotplot(kk, showCategory = top_n) + ggtitle(paste("KEGG ORA:", dir)), width = 9, height = 6, dpi = 150) } else empty(fk)
}
if (!run_kegg) kegg_src <- "KEGG skipped (--run_kegg false, the default; KEGG licence for commercial use is unresolved)"
writeLines(kegg_src, "ora_kegg_source.txt"); message(kegg_src)
