#!/usr/bin/env Rscript
# Reference-gene check on the 19 diagnostic samples (featureCounts). Run from the project root with env/prefix/r_dge/bin/Rscript. The 15 validation samples are not used.
# Start list: Eisenberg & Levanon 2013 housekeeping genes (results/13_reference_genes/source/HK_genes.txt, md5 3c29a7d8976857b6242b317bf58ba25b).
# Order: map -> keep genes tested on both routes -> remove rule-1 genes (either route) -> baseMean >= 100 on both routes (PROPOSED) -> VST stats, log2FC with 90% CI, study coefficient -> equivalence.
suppressPackageStartupMessages({ library(DESeq2) })
out <- "results/13_reference_genes"; BM_MIN <- 100; CI_LEVEL <- 0.90; EQ <- c(0.5, 0.3)
strip <- function(x) sub("\\..*$", "", x)
hk <- read.table(file.path(out, "source/HK_genes.txt"), sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, strip.white = TRUE, col.names = c("symbol", "refseq"), fill = TRUE)
gs <- read.delim("ref/gene_id_to_symbol.tsv", stringsAsFactors = FALSE); gs$id <- strip(gs$gene_id)
resc <- read.delim(file.path(out, "hk_refseq_rescue_map.tsv"), stringsAsFactors = FALSE)
fc <- read.delim("results/06_dge/deseq2_results_all.tsv", stringsAsFactors = FALSE); sa <- read.delim("results/06_dge/deseq2_salmon_results_all.tsv", stringsAsFactors = FALSE)
fc$id <- strip(fc$gene_id); sa$id <- strip(sa$gene_id)
tested <- intersect(fc$id[!is.na(fc$pvalue)], sa$id[!is.na(sa$pvalue)])
ct <- read.delim("results/11_candidate_selection/candidate_table.tsv", stringsAsFactors = FALSE); rule1_ids <- strip(ct$gene_id[ct$rule1_pass])   # rule 1 on either route

# ---- 2. mapping: (a) exact symbol -> GENCODE v50 gene_name; (b) RefSeq accession -> Ensembl (org.Hs.eg.db) for symbols not found in (a). Genes with >1 Ensembl ID for one symbol are ambiguous and dropped.
a <- unique(merge(hk, gs[, c("gene_name", "id")], by.x = "symbol", by.y = "gene_name")[, c("symbol", "id")]); a$method <- "symbol_GENCODE_v50"
b <- unique(data.frame(symbol = resc$symbol, id = resc$ENSEMBL, method = "RefSeq_orgHsEgDb", stringsAsFactors = FALSE)); b <- b[!(b$symbol %in% a$symbol), ]
mp <- rbind(a, b); amb <- names(which(table(unique(mp[, c("symbol", "id")])$symbol) > 1)); mp_amb <- mp[mp$symbol %in% amb, ]; mp <- mp[!(mp$symbol %in% amb), ]
mp <- mp[!duplicated(mp$id), ]
steps <- data.frame(step = character(), genes = integer(), note = character(), stringsAsFactors = FALSE)
add <- function(s, n, note = "") steps[nrow(steps) + 1, ] <<- list(s, as.integer(n), note)
add("HK list (Eisenberg & Levanon 2013), unique symbols", length(unique(hk$symbol)))
add("symbols matched to GENCODE v50 gene_name (exact)", length(unique(a$symbol)), "")
add("symbols added by RefSeq accession (org.Hs.eg.db), not found by symbol", length(unique(b$symbol)), "")
add("ambiguous symbols (>1 Ensembl ID), dropped", length(amb))
add("symbols with no Ensembl ID by either method (unmatched)", length(setdiff(hk$symbol, c(a$symbol, b$symbol))))
add("mapped Ensembl IDs (unique, unversioned)", nrow(mp))
mp$tested_both <- mp$id %in% tested; add("  of which tested on both routes", sum(mp$tested_both)); k <- mp[mp$tested_both, ]
k$rule1_either <- k$id %in% rule1_ids; add("removed: pass rule 1 on either route", sum(k$rule1_either)); k <- k[!k$rule1_either, ]; add("after removing rule-1 genes", nrow(k))
k$bm_fc <- fc$baseMean[match(k$id, fc$id)]; k$bm_sal <- sa$baseMean[match(k$id, sa$id)]
k <- k[k$bm_fc >= BM_MIN & k$bm_sal >= BM_MIN, ]; add("after baseMean >= 100 on both routes (PROPOSED)", nrow(k))

# ---- 19-sample fit (same call as before) for study coefficient and unshrunk group LFC; VST over all 19 (blind), as arm 1 of the LASSO trial
cn <- round(as.matrix(read.delim("results/04_counts/counts_featurecounts.tsv", row.names = 1, check.names = FALSE))); rownames(cn) <- rownames(cn)
sm <- read.delim("results/sample_map.tsv", stringsAsFactors = FALSE); rownames(sm) <- sm$srr; ss <- sm[colnames(cn), ]; ss$study <- factor(ss$study); ss$group <- relevel(factor(ss$group), ref = "normal")
dds <- DESeqDataSetFromMatrix(cn, colData = ss, design = ~ study + group); dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]; dds <- DESeq(dds, quiet = TRUE)
print(resultsNames(dds)); stopifnot(all(c("group_tumour_vs_normal", "study_SRP301216_vs_SRP219837") %in% resultsNames(dds)))
rg <- results(dds, name = "group_tumour_vs_normal"); rs <- results(dds, name = "study_SRP301216_vs_SRP219837")
rownames(rg) <- strip(rownames(rg)); rownames(rs) <- strip(rownames(rs))
vst <- assay(varianceStabilizingTransformation(dds, blind = TRUE)); rownames(vst) <- strip(rownames(vst)); nc <- counts(dds, normalized = TRUE); rownames(nc) <- strip(rownames(nc))
chk <- max(abs(rg[k$id, "log2FoldChange"] - fc$log2FoldChange[match(k$id, fc$id)]), na.rm = TRUE); message("max abs diff group LFC refit vs saved table: ", signif(chk, 3)); stopifnot(chk < 1e-4)
z <- qnorm(1 - (1 - CI_LEVEL) / 2)
cv <- function(x) sd(x) / mean(x); stA <- ss$study == "SRP219837"; stB <- !stA
k$symbol <- k$symbol; k$gene_biotype <- gs$gene_type[match(k$id, gs$id)]
k$vst_mean <- rowMeans(vst[k$id, ]); k$vst_cv <- apply(vst[k$id, ], 1, cv)
k$vst_mean_A <- rowMeans(vst[k$id, stA]); k$vst_cv_A <- apply(vst[k$id, stA], 1, cv); k$vst_mean_B <- rowMeans(vst[k$id, stB]); k$vst_cv_B <- apply(vst[k$id, stB], 1, cv)
k$normcount_cv <- apply(nc[k$id, ], 1, cv)                                   # extra: CV on the count scale (the VST CV is on a log scale and is compressed)
k$log2FC <- rg[k$id, "log2FoldChange"]; k$lfcSE <- rg[k$id, "lfcSE"]; k$ci90_low <- k$log2FC - z * k$lfcSE; k$ci90_high <- k$log2FC + z * k$lfcSE; k$padj_group <- rg[k$id, "padj"]
k$study_coef_B_vs_A <- rs[k$id, "log2FoldChange"]; k$study_coef_se <- rs[k$id, "lfcSE"]
k$log2FC_salmon <- sa$log2FoldChange[match(k$id, sa$id)]; k$lfcSE_salmon <- sa$lfcSE[match(k$id, sa$id)]
for (e in EQ) k[[paste0("equivalent_", e)]] <- k$ci90_low > -e & k$ci90_high < e
add("equivalent tumour vs normal, 90% CI inside +/-0.5 log2", sum(k$equivalent_0.5)); add("equivalent, 90% CI inside +/-0.3 log2", sum(k$equivalent_0.3))
k <- k[order(!k$equivalent_0.5, k$vst_cv), ]; k$rank <- seq_len(nrow(k))
k <- k[, c("rank", "id", "symbol", "gene_biotype", "method", "bm_fc", "bm_sal", "vst_mean", "vst_cv", "vst_mean_A", "vst_cv_A", "vst_mean_B", "vst_cv_B", "normcount_cv", "log2FC", "lfcSE", "ci90_low", "ci90_high", "padj_group", "study_coef_B_vs_A", "study_coef_se", "log2FC_salmon", "lfcSE_salmon", "equivalent_0.5", "equivalent_0.3")]
write.table(k, file.path(out, "reference_gene_candidates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(steps, file.path(out, "step_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mp_amb, file.path(out, "hk_ambiguous_symbols.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(symbol = setdiff(hk$symbol, c(a$symbol, b$symbol))), file.path(out, "hk_unmatched_symbols.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(steps, right = FALSE); top <- head(k, 10); print(top[, c("rank", "symbol", "vst_mean", "vst_cv", "log2FC", "ci90_low", "ci90_high", "study_coef_B_vs_A", "equivalent_0.5", "equivalent_0.3")], digits = 3)
writeLines(c(paste("R", getRversion()), paste("DESeq2", packageVersion("DESeq2")), paste("run", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), file.path(out, "tool_versions.txt"))
