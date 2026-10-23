#!/usr/bin/env Rscript
# NEW. Compare with Hosseini et al. 2023 expectations AFTER the analysis. Never feeds back into cutoffs.
# usage: check_vs_paper.R <deseq2_results_all.tsv> <ora_dir> <gsea_dir> <padj>
a <- commandArgs(trailingOnly = TRUE); pc <- as.numeric(a[4])
d <- read.delim(a[1], stringsAsFactors = FALSE)
down <- c("GUCA2A", "CLCA4", "CLCA1", "AQP8", "SLC26A3", "APOA4", "PYY", "GUCA2B", "ZG16", "CD177")
up   <- c("COL3A1", "COL1A1", "SPARC", "TGFBI", "CXCL8", "KRT6A", "KRT6B", "KRT17")
g <- data.frame(gene = c(down, up), expected = c(rep("down", length(down)), rep("up", length(up))), stringsAsFactors = FALSE)
r <- do.call(rbind, lapply(g$gene, function(x) { z <- d[!is.na(d$gene_name) & d$gene_name == x, ]
  if (!nrow(z)) data.frame(log2FC = NA_real_, log2FC_shrunk = NA_real_, padj = NA_real_)
  else { z <- z[which.max(z$baseMean), ]; data.frame(log2FC = z$log2FoldChange, log2FC_shrunk = z$log2FoldChange_shrunk, padj = z$padj) } }))
g <- cbind(g, r)
g$observed <- ifelse(is.na(g$log2FC), "not tested", ifelse(g$log2FC > 0, "up", "down"))
g$significant <- !is.na(g$padj) & g$padj < pc
g$match_direction <- g$observed == g$expected
g$match_direction_and_significant <- g$match_direction & g$significant
write.table(g, "check_vs_paper_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

norm <- function(x) { x <- tolower(x); x <- sub(" - homo sapiens \\(human\\)$", "", x); x <- sub("^kegg[_ ]", "", x); x <- gsub("_", " ", x)
  x <- gsub("-", "", x); x <- sub(" pathway$", "", x); trimws(x) }   # e.g. KEGG_IL17_SIGNALING_PATHWAY and "IL-17 signaling pathway" -> "il17 signaling"
pw <- data.frame(pathway = norm(c("mineral absorption", "fat digestion and absorption", "ppar signaling", "pancreatic secretion", "nitrogen metabolism", "bile secretion",
                             "cell cycle", "spliceosome", "rna transport", "dna replication", "mrna surveillance", "il-17 signaling")),
                 expected = c(rep("down", 6), rep("up", 6)), stringsAsFactors = FALSE)
rd <- function(f) if (file.exists(f) && file.info(f)$size > 0) read.delim(f, stringsAsFactors = FALSE) else NULL
find_hit <- function(df, col, p) { if (is.null(df) || !nrow(df)) return(which(FALSE)); which(grepl(p, norm(df[[col]]), fixed = TRUE)) }
ora_up <- rd(file.path(a[2], "ora_up_KEGG.tsv")); ora_dn <- rd(file.path(a[2], "ora_down_KEGG.tsv")); gs <- rd(file.path(a[3], "gsea_KEGG.tsv"))
pw$ora_up_padj <- pw$ora_down_padj <- pw$gsea_NES <- pw$gsea_padj <- NA_real_
for (i in seq_len(nrow(pw))) { p <- pw$pathway[i]
  h <- find_hit(ora_up, "Description", p); if (!length(h)) h <- find_hit(ora_up, "ID", p); if (length(h)) pw$ora_up_padj[i] <- min(ora_up$p.adjust[h])
  h <- find_hit(ora_dn, "Description", p); if (!length(h)) h <- find_hit(ora_dn, "ID", p); if (length(h)) pw$ora_down_padj[i] <- min(ora_dn$p.adjust[h])
  h <- find_hit(gs, "pathway", p); if (length(h)) { h <- h[which.min(gs$padj[h])]; pw$gsea_NES[i] <- gs$NES[h]; pw$gsea_padj[i] <- gs$padj[h] } }
sigdir <- function(i) { dirs <- c(); if (!is.na(pw$ora_up_padj[i]) && pw$ora_up_padj[i] < pc) dirs <- c(dirs, "up")
  if (!is.na(pw$ora_down_padj[i]) && pw$ora_down_padj[i] < pc) dirs <- c(dirs, "down")
  if (!is.na(pw$gsea_padj[i]) && pw$gsea_padj[i] < pc) dirs <- c(dirs, ifelse(pw$gsea_NES[i] > 0, "up", "down")); dirs }
pw$observed_significant_directions <- vapply(seq_len(nrow(pw)), function(i) { s <- unique(sigdir(i)); if (!length(s)) "none significant" else paste(s, collapse = "+") }, "")
pw$match <- mapply(function(e, o) o == e, pw$expected, pw$observed_significant_directions)
write.table(pw, "check_vs_paper_pathways.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
txt <- c(sprintf("Genes: direction matches %d of %d; direction matches AND padj < %g: %d of %d; not tested: %d", sum(g$match_direction), nrow(g), pc, sum(g$match_direction_and_significant), nrow(g), sum(g$observed == "not tested")),
         sprintf("  expected down: %d of %d match direction; expected up: %d of %d match direction", sum(g$match_direction[g$expected == "down"]), sum(g$expected == "down"), sum(g$match_direction[g$expected == "up"]), sum(g$expected == "up")),
         sprintf("Pathways (KEGG; ORA or GSEA significant at q < %g in the expected direction): %d of %d", pc, sum(pw$match), nrow(pw)),
         "A mismatch is information, not a reason to change cutoffs (paper used 4 SRA studies + TCGA; overlap expected, not identity).")
writeLines(txt, "check_vs_paper_summary.txt"); cat(txt, sep = "\n")
