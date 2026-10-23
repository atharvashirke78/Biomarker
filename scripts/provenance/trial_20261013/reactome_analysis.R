#!/usr/bin/env Rscript
# Reactome via MSigDB C2:CP:REACTOME (msigdbr), NOT the official Reactome release. Run from the project root:
#   R_USER_CACHE_DIR=$PWD/env/msigdbr_cache env/prefix/r_pathway/bin/Rscript results/12_reactome/reactome_analysis.R
# Reads only: results/11_candidate_selection/candidate_table.tsv, results/06_dge/deseq2_results_all.tsv, deseq2_salmon_results_all.tsv. Writes only here.
# Gene IDs: unversioned Ensembl gene ID on both sides (ours: version/_PAR_Y suffix stripped; MSigDB: ensembl_gene). No symbol matching.
suppressPackageStartupMessages({ library(msigdbr); library(clusterProfiler); library(fgsea); library(BiocParallel) })
out <- "results/12_reactome"
MIN_GS <- 10; MAX_GS <- 500; Q <- 0.05; NPERM <- 10000; SEED <- 42      # same as conf/thresholds.config:43-47 (path_min_gs, path_max_gs, path_q, gsea_nperm, gsea_seed)
strip <- function(x) sub("\\..*$", "", x)

rs <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
stopifnot(nrow(rs) > 0)
t2g <- unique(data.frame(term = rs$gs_name, gene = rs$ensembl_gene, stringsAsFactors = FALSE)); t2g <- t2g[!is.na(t2g$gene) & t2g$gene != "", ]
db_version <- unique(rs$db_version); n_sets <- length(unique(t2g$term))

fc <- read.delim("results/06_dge/deseq2_results_all.tsv", stringsAsFactors = FALSE); sa <- read.delim("results/06_dge/deseq2_salmon_results_all.tsv", stringsAsFactors = FALSE)
fc$id <- strip(fc$gene_id); sa$id <- strip(sa$gene_id)
dup <- c(featureCounts = sum(duplicated(fc$id)), salmon = sum(duplicated(sa$id))); print(dup)
fc <- fc[!duplicated(fc$id), ]; sa <- sa[!duplicated(sa$id), ]
bg <- intersect(fc$id[!is.na(fc$pvalue)], sa$id[!is.na(sa$pvalue)])            # tested (non-NA p-value) on both routes
ct <- read.delim("results/11_candidate_selection/candidate_table.tsv", stringsAsFactors = FALSE); ct$id <- strip(ct$gene_id)
cand <- ct[ct$all_rules_pass, ]
up <- unique(cand$id[cand$direction_fc == "up"]); down <- unique(cand$id[cand$direction_fc == "down"]); both <- unique(cand$id)
inR <- unique(t2g$gene)
ms <- data.frame(item = c("MSigDB release (db_version)", "msigdbr version", "gene sets in C2:CP:REACTOME", "unique Ensembl IDs in the subset",
    "featureCounts genes tested (non-NA p)", "  of which in Reactome subset", "background: tested on both routes", "  of which in Reactome subset",
    "candidates after rule 3b", "  of which in background", "  of which in Reactome subset", "  up / down (featureCounts direction)", "duplicate unversioned IDs removed (fc / salmon)"),
  value = c(db_version, as.character(packageVersion("msigdbr")), n_sets, length(inR),
    sum(!is.na(fc$pvalue)), sum(fc$id[!is.na(fc$pvalue)] %in% inR), length(bg), sum(bg %in% inR),
    length(both), sum(both %in% bg), sum(both %in% inR), paste(length(up), length(down), sep = " / "), paste(dup, collapse = " / ")))
print(ms); write.table(ms, file.path(out, "mapping_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
id2sym <- setNames(fc$gene_name, fc$id)

ora <- function(g, name) {
  g <- intersect(g, bg)
  e <- enricher(g, universe = bg, TERM2GENE = t2g, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS, maxGSSize = MAX_GS)
  r <- if (is.null(e)) data.frame() else as.data.frame(e)
  if (nrow(r)) { r$significant_BH_lt_0.05 <- r$p.adjust < Q
    r$geneSymbols <- vapply(strsplit(r$geneID, "/"), function(x) paste(id2sym[x], collapse = ";"), "") }
  r$input_genes_in_background <- if (nrow(r)) length(g) else integer(0)
  write.table(r, file.path(out, paste0("ora_", name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message(name, ": ", length(g), " genes, ", nrow(r), " sets tested, ", sum(r$significant_BH_lt_0.05), " with BH < ", Q); r }
r_up <- ora(up, "up"); r_down <- ora(down, "down"); r_both <- ora(both, "combined")

d <- fc[!is.na(fc$stat), ]; ranks <- sort(setNames(d$stat, d$id), decreasing = TRUE)      # Wald statistic, featureCounts, one value per unversioned ID
sets <- split(t2g$gene, t2g$term)
set.seed(SEED)
fg <- fgseaMultilevel(sets, ranks, minSize = MIN_GS, maxSize = MAX_GS, nPermSimple = NPERM, BPPARAM = SerialParam())
fg <- fg[order(fg$padj), ]; fg$leadingEdge <- vapply(fg$leadingEdge, paste, character(1), collapse = ";"); fg$significant_padj_lt_0.05 <- !is.na(fg$padj) & fg$padj < Q
write.table(as.data.frame(fg), file.path(out, "gsea_reactome.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("GSEA: ", length(ranks), " ranked genes, ", nrow(fg), " sets tested, ", sum(fg$significant_padj_lt_0.05), " with padj < ", Q)
write.table(data.frame(analysis = c("ORA up", "ORA down", "ORA combined", "GSEA"), genes_in = c(length(intersect(up, bg)), length(intersect(down, bg)), length(intersect(both, bg)), length(ranks)),
  sets_tested = c(nrow(r_up), nrow(r_down), nrow(r_both), nrow(fg)), sets_significant = c(sum(r_up$significant_BH_lt_0.05), sum(r_down$significant_BH_lt_0.05), sum(r_both$significant_BH_lt_0.05), sum(fg$significant_padj_lt_0.05))),
  file.path(out, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("R", getRversion()), paste("msigdbr", packageVersion("msigdbr")), paste("MSigDB db_version", db_version), paste("clusterProfiler", packageVersion("clusterProfiler")), paste("fgsea", packageVersion("fgsea")),
  "source: Reactome via MSigDB C2:CP:REACTOME (msigdbr), NOT the official Reactome release", paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
