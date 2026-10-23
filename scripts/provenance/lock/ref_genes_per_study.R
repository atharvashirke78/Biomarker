# DRAFT v2: tumour vs normal log2FC (unshrunk DESeq2 Wald) with 90% CI for the 10 frozen reference genes, WITHIN each study, 19 development samples only.
# Model ~ group within one study (reference normal); pre-filter >= 10 counts in >= n samples, n = smallest group of THAT study (same 'auto' rule); CI = log2FC +/- qnorm(0.95) * lfcSE (as in reference_gene_check.R).
suppressPackageStartupMessages(library(DESeq2))
cn <- round(as.matrix(read.delim("results/04_counts/counts_featurecounts.tsv", row.names = 1, check.names = FALSE))); stopifnot(ncol(cn) == 19)
sm <- read.delim("results/sample_map.tsv", stringsAsFactors = FALSE); rownames(sm) <- sm$srr; ss <- sm[colnames(cn), ]
genes <- read.delim("results_rc1prep_19/lock_proposal/reference_gene_list_top10.txt", header = FALSE, col.names = c("symbol", "id"), stringsAsFactors = FALSE)
MIN_COUNT <- 10; z <- qnorm(0.95); rows <- list()
for (st in sort(unique(ss$study))) {
  cols <- rownames(ss)[ss$study == st]; cd <- ss[cols, ]; cd$group <- relevel(factor(cd$group), ref = "normal"); n_min <- min(table(cd$group))
  d <- DESeqDataSetFromMatrix(cn[, cols], colData = cd, design = ~ group); d <- d[rowSums(counts(d) >= MIN_COUNT) >= n_min, ]; d <- DESeq(d, quiet = TRUE)
  r <- results(d, name = "group_tumour_vs_normal"); rownames(r) <- sub("\\..*$", "", rownames(r))
  cat(st, ": samples", length(cols), "(", paste(names(table(cd$group)), table(cd$group), collapse = ", "), ") pre-filter n =", n_min, ", genes tested", sum(!is.na(r$pvalue)), "\n")
  m <- match(genes$id, rownames(r))
  rows[[st]] <- data.frame(study = st, symbol = genes$symbol, id = genes$id, n_normal = sum(cd$group == "normal"), n_tumour = sum(cd$group == "tumour"), baseMean = r$baseMean[m], log2FC = r$log2FoldChange[m], lfcSE = r$lfcSE[m],
    ci90_low = r$log2FoldChange[m] - z * r$lfcSE[m], ci90_high = r$log2FoldChange[m] + z * r$lfcSE[m], stringsAsFactors = FALSE)
}
T <- do.call(rbind, rows); T$inside_0.5 <- T$ci90_low > -0.5 & T$ci90_high < 0.5; T$inside_0.3 <- T$ci90_low > -0.3 & T$ci90_high < 0.3
write.table(T, "results_rc1prep_19/lock_proposal_v2/ref_genes_per_study.tsv", sep = "\t", quote = FALSE, row.names = FALSE); print(T, digits = 3)
