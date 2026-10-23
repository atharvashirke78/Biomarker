#!/usr/bin/env Rscript
# Step 10: significant set redefined as padj < 0.05 & |apeglm-shrunk log2FC| >= 1 (pipeline uses unshrunk). Reads results only; nothing in the pipeline or earlier results is changed.
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
B <- "/workspace/Biomarker_pipeline"; R <- file.path(B, "results"); V <- file.path(R, "09_verification"); O <- file.path(V, "step10_shrunk_cutoff"); X <- file.path(V, "input/MOESM2_sheets")
PADJ <- 0.05; LFC <- 1; `%+%` <- function(a, b) paste0(a, b)
sink(file.path(O, "step10_log.txt"), split = TRUE)
rd0 <- function(f) read.delim(file.path(R, "06_dge", f), stringsAsFactors = FALSE)
res <- list(featureCounts = rd0("deseq2_results_all.tsv"), Salmon = rd0("deseq2_salmon_results_all.tsv"))
sig_u <- function(o) !is.na(o$padj) & o$padj < PADJ & abs(o$log2FoldChange) >= LFC
sig_s <- function(o) !is.na(o$padj) & o$padj < PADJ & abs(o$log2FoldChange_shrunk) >= LFC
# ---------------- 1. sets
s1 <- c("== 1. Significant set: padj < 0.05 & |log2FC| >= 1, unshrunk (pipeline) vs apeglm-shrunk", "")
S1 <- data.frame()
for (rt in names(res)) { o <- res[[rt]]; u <- o$gene_id[sig_u(o)]; s <- o$gene_id[sig_s(o)]; ov <- intersect(u, s)
  lost <- setdiff(u, s); gained <- setdiff(s, u)
  su <- sig_s(o)
  S1 <- rbind(S1, data.frame(route = rt, tested = sum(!is.na(o$padj)), unshrunk_n = length(u), unshrunk_up = sum(sig_u(o) & o$log2FoldChange > 0), unshrunk_down = sum(sig_u(o) & o$log2FoldChange < 0),
    shrunk_n = length(s), shrunk_up = sum(su & o$log2FoldChange_shrunk > 0), shrunk_down = sum(su & o$log2FoldChange_shrunk < 0), overlap = length(ov), lost = length(lost), gained = length(gained), jaccard = round(length(ov) / length(union(u, s)), 4)))
  write.table(o[o$gene_id %in% lost, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "log2FoldChange_shrunk", "lfcSE", "padj")], file.path(O, sprintf("genes_lost_%s.tsv", rt)), sep = "\t", quote = FALSE, row.names = FALSE)
  l <- o[o$gene_id %in% lost, ]
  s1 <- c(s1, sprintf("%s: shrunk set %d genes (up %d, down %d) vs unshrunk %d (up %d, down %d); overlap %d; Jaccard %.3f; lost %d (%.1f%% of unshrunk); gained %d",
     rt, length(s), sum(su & o$log2FoldChange_shrunk > 0), sum(su & o$log2FoldChange_shrunk < 0), length(u), sum(sig_u(o) & o$log2FoldChange > 0), sum(sig_u(o) & o$log2FoldChange < 0), length(ov), length(ov) / length(union(u, s)), length(lost), 100 * length(lost) / length(u), length(gained)),
     sprintf("   lost genes: median baseMean %.1f (kept genes %.1f); median unshrunk |lfc| %.2f -> shrunk %.2f; median lfcSE %.2f (kept %.2f); direction flips (shrunk sign != unshrunk sign): %d",
        median(l$baseMean), median(o$baseMean[o$gene_id %in% ov]), median(abs(l$log2FoldChange)), median(abs(l$log2FoldChange_shrunk)), median(l$lfcSE), median(o$lfcSE[o$gene_id %in% ov]), sum(sign(l$log2FoldChange) != sign(l$log2FoldChange_shrunk))), "")
}
write.table(S1, file.path(O, "sets_unshrunk_vs_shrunk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# route concordance of the shrunk sets
sf <- res$featureCounts$gene_id[sig_s(res$featureCounts)]; ss_ <- res$Salmon$gene_id[sig_s(res$Salmon)]; uf <- res$featureCounts$gene_id[sig_u(res$featureCounts)]; us <- res$Salmon$gene_id[sig_u(res$Salmon)]
s1 <- c(s1, sprintf("featureCounts vs Salmon overlap (gene_id): unshrunk %d common of %d / %d (Jaccard %.3f); shrunk %d common of %d / %d (Jaccard %.3f)", length(intersect(uf, us)), length(uf), length(us), length(intersect(uf, us)) / length(union(uf, us)),
     length(intersect(sf, ss_)), length(sf), length(ss_), length(intersect(sf, ss_)) / length(union(sf, ss_))), "")
writeLines(s1, file.path(O, "part1_sets_summary.txt")); cat(s1, sep = "\n")
# ---------------- 2. Table S3 and 18 genes
S3 <- read.delim(file.path(V, "input/hosseini_2023_TableS3_genes.tsv"), stringsAsFactors = FALSE); S3$gene <- trimws(S3$gene); S3 <- S3[!duplicated(S3$gene), ]
s2 <- c("== 2. Hosseini Table S3 agreement (same method as step5; universe = tested genes with symbol, one row per symbol = highest baseMean)", "")
for (rt in names(res)) {
  F <- res[[rt]]; Uo <- F[!is.na(F$gene_name) & F$gene_name != "" & !is.na(F$padj), ]; Uo <- Uo[order(-Uo$baseMean), ]; Uo <- Uo[!duplicated(Uo$gene_name), ]; N <- nrow(Uo)
  m0 <- merge(Uo[, c("gene_name", "log2FoldChange", "log2FoldChange_shrunk", "padj")], S3, by.x = "gene_name", by.y = "gene"); m0$ours <- ifelse(m0$log2FoldChange > 0, "up", "down"); m0$agree <- m0$ours == m0$direction_in_cancer
  for (def in c("unshrunk", "shrunk")) {
    ok <- if (def == "unshrunk") sig_u(Uo) else sig_s(Uo); lf <- if (def == "unshrunk") Uo$log2FoldChange else Uo$log2FoldChange_shrunk
    su <- Uo$gene_name[ok & lf > 0]; sd <- Uo$gene_name[ok & lf < 0]
    sigset <- Uo$gene_name[ok]; m <- m0; m$sig <- m$gene_name %in% sigset
    pu <- m$gene_name[m$direction_in_cancer == "up"]; pd <- m$gene_name[m$direction_in_cancer == "down"]
    hyp <- function(a, b) { k <- length(intersect(a, b)); K <- length(b); n <- length(a); c(ours = n, paper = K, overlap = k, expected = round(n * K / N, 1), fold = round(k / (n * K / N), 2), p = signif(phyper(k - 1, K, N - K, n, lower.tail = FALSE), 3)) }
    H <- rbind(`ours UP vs S3 up` = hyp(su, pu), `ours DOWN vs S3 down` = hyp(sd, pd), `ours UP vs S3 down (discordant)` = hyp(su, pd), `ours DOWN vs S3 up (discordant)` = hyp(sd, pu))
    dsum <- function(d) { x <- m[m$direction_in_cancer == d, ]; y <- x[x$sig, ]; sprintf("   S3 %-4s: %d tested; sign agrees %d (%.1f%%); among our significant: %d, agree %d (%.1f%%)", d, nrow(x), sum(x$agree), 100 * mean(x$agree), nrow(y), sum(y$agree), 100 * mean(y$agree)) }
    s2 <- c(s2, sprintf("-- %s, %s significant set (%d genes with symbol)", rt, def, length(sigset)), dsum("up"), dsum("down"), "   hypergeometric (N = " %+% N %+% "):", paste0("   ", capture.output(print(H))), "")
    write.table(cbind(test = rownames(H), H), file.path(O, sprintf("S3_hypergeometric_%s_%s.tsv", rt, def)), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
g18 <- c("GUCA2A","CLCA4","CLCA1","AQP8","SLC26A3","APOA4","PYY","GUCA2B","ZG16","CD177","COL3A1","COL1A1","SPARC","TGFBI","CXCL8","KRT6A","KRT6B","KRT17"); e18 <- c(rep("down", 10), rep("up", 8))
for (rt in names(res)) {
  F <- res[[rt]]; t <- do.call(rbind, lapply(seq_along(g18), function(i) { z <- F[!is.na(F$gene_name) & F$gene_name == g18[i], ]
    if (!nrow(z)) return(data.frame(gene = g18[i], expected = e18[i], lfc = NA, lfc_shrunk = NA, padj = NA, sig_unshrunk = NA, sig_shrunk = NA, direction = NA, dir_ok = NA))
    z <- z[which.max(z$baseMean), ]; d <- ifelse(z$log2FoldChange > 0, "up", "down")
    data.frame(gene = g18[i], expected = e18[i], lfc = z$log2FoldChange, lfc_shrunk = z$log2FoldChange_shrunk, padj = z$padj, sig_unshrunk = sig_u(z), sig_shrunk = sig_s(z), direction = d, dir_ok = d == e18[i]) }))
  write.table(t, file.path(O, sprintf("hosseini18_%s.tsv", rt)), sep = "\t", quote = FALSE, row.names = FALSE)
  tt <- t; tt$lfc <- round(tt$lfc, 2); tt$lfc_shrunk <- round(tt$lfc_shrunk, 2); tt$padj <- signif(tt$padj, 2)
  s2 <- c(s2, sprintf("-- 18 genes, %s (padj NA = Cook's outlier flag or not tested; sig_* needs padj < 0.05 and |lfc| >= 1): direction ok %d/%d testable; significant unshrunk %d, shrunk %d", rt, sum(t$dir_ok, na.rm = TRUE), sum(!is.na(t$dir_ok)), sum(t$sig_unshrunk, na.rm = TRUE), sum(t$sig_shrunk, na.rm = TRUE)), capture.output(print(tt, row.names = FALSE)), "")
}
writeLines(s2, file.path(O, "part2_S3_18genes_summary.txt")); cat(s2, sep = "\n")
# ---------------- 3. ORA (same settings as pipeline ora.R: universe = tested genes with symbol and pvalue, set size 10-500, BH q 0.05) for the shrunk set, featureCounts route
d <- res$featureCounts; d <- d[!is.na(d$gene_name) & d$gene_name != "" & !is.na(d$pvalue), ]
sg <- !is.na(d$padj) & d$padj < PADJ & abs(d$log2FoldChange_shrunk) >= LFC
sets <- list(up = unique(d$gene_name[sg & d$log2FoldChange_shrunk > 0]), down = unique(d$gene_name[sg & d$log2FoldChange_shrunk < 0]))
map <- bitr(unique(d$gene_name), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db); uni_e <- unique(map$ENTREZID)
rd <- function(f) read.delim(file.path(X, f), skip = 1, stringsAsFactors = FALSE, quote = "")
stripgo <- function(x) trimws(sub("\\s*\\(GO:[0-9]+\\)\\s*$", "", x)); norm <- function(x) tolower(gsub("[^a-z0-9]+", " ", tolower(trimws(x))))
kegg_ok <- TRUE; rows <- list()
for (dir in c("up", "down")) {
  g <- sets[[dir]]; ge <- unique(map$ENTREZID[map$SYMBOL %in% g]); message(dir, ": ", length(g), " symbols, ", length(ge), " Entrez ids")
  go <- enrichGO(ge, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = uni_e, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, readable = TRUE)
  kk <- tryCatch(enrichKEGG(ge, organism = "hsa", universe = uni_e, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500), error = function(e) { message("enrichKEGG failed: ", conditionMessage(e)); NULL })
  if (is.null(kk)) kegg_ok <- FALSE else kk <- setReadable(kk, org.Hs.eg.db, keyType = "ENTREZID")
  goT <- as.data.frame(go); write.table(goT, file.path(O, sprintf("ora_shrunk_%s_GO_BP.tsv", dir)), sep = "\t", quote = FALSE, row.names = FALSE)
  kT <- if (is.null(kk)) data.frame() else as.data.frame(kk); write.table(kT, file.path(O, sprintf("ora_shrunk_%s_KEGG.tsv", dir)), sep = "\t", quote = FALSE, row.names = FALSE)
  for (typ in c("KEGG", "GO BP")) {
    paper <- if (typ == "KEGG") { p <- rd(if (dir == "up") "Table_S6_Upregulated_.tsv" else "Table_S6_Downregulated_.tsv"); p[p$Adjusted.P.value < 0.05, "Term"] } else { p <- rd(if (dir == "up") "Table_S5_Biological_process_.tsv" else "Table_S4_Biological_process_.tsv"); stripgo(p[p$Adjusted.P.value < 0.05, "Term"]) }
    old <- read.delim(file.path(R, "07_pathway/ora", sprintf("ora_%s_%s.tsv", dir, sub(" ", "_", typ))), stringsAsFactors = FALSE); old <- old[old$p.adjust < 0.05, ]
    new <- if (typ == "KEGG") kT else goT; new <- if (nrow(new)) new[new$p.adjust < 0.05, ] else new
    pn <- unique(norm(paper)); on <- unique(norm(old$Description)); nn <- if (nrow(new)) unique(norm(new$Description)) else character()
    rows[[length(rows) + 1]] <- data.frame(set = paste(typ, dir), paper = length(pn), unshrunk_ours = length(on), unshrunk_shared = length(intersect(on, pn)), shrunk_ours = length(nn), shrunk_shared = length(intersect(nn, pn)),
       shrunk_paper_only = length(setdiff(pn, nn)), shrunk_ours_only = length(setdiff(nn, pn)), shared_lost_vs_unshrunk = length(setdiff(intersect(on, pn), nn)), shared_gained = length(setdiff(intersect(nn, pn), on)))
  }
}
T3 <- do.call(rbind, rows); write.table(T3, file.path(O, "ora_overlap_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
s3 <- c("== 3. ORA terms shared with the paper (padj < 0.05 both; names matched after lower-casing and removing punctuation), featureCounts route",
  sprintf("shrunk gene sets: up %d, down %d symbols (unshrunk: up %d, down %d); KEGG source: %s", length(sets$up), length(sets$down), length(unique(d$gene_name[!is.na(d$padj) & d$padj < PADJ & d$log2FoldChange >= LFC])), length(unique(d$gene_name[!is.na(d$padj) & d$padj < PADJ & d$log2FoldChange <= -LFC])),
          if (kegg_ok) "KEGG REST (live release at run time; the unshrunk ORA used the release of 2026-10-09, so KEGG differences can come from the release too)" else "KEGG REST FAILED, KEGG columns empty"),
  "unshrunk = existing results/07_pathway/ora files, not rerun.", "", capture.output(print(T3, row.names = FALSE)))
writeLines(s3, file.path(O, "part3_ora_summary.txt")); cat(s3, sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(O, "step10_sessionInfo.txt")); sink()
