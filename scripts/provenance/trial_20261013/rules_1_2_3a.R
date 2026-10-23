# rules 1, 2, 3a only (same definitions as results/11_candidate_selection/candidate_selection.R lines 30, 49-62); reads results_rc1prep_19/deseq2*_results_all.tsv; no 3b, no LOO
fc <- read.delim("results_rc1prep_19/deseq2_results_all.tsv", stringsAsFactors = FALSE); sa <- read.delim("results_rc1prep_19/deseq2_salmon_results_all.tsv", stringsAsFactors = FALSE)
r1 <- function(p, l) !is.na(p) & p < 0.05 & !is.na(l) & abs(l) >= 1
fc$r1 <- r1(fc$padj, fc$log2FoldChange_shrunk); sa$r1 <- r1(sa$padj, sa$log2FoldChange_shrunk)
ids <- union(fc$gene_id[fc$r1], sa$gene_id[sa$r1]); f <- fc[match(ids, fc$gene_id), ]; s <- sa[match(ids, sa$gene_id), ]
t <- data.frame(gene_id = ids, symbol = ifelse(!is.na(f$gene_name), f$gene_name, s$gene_name), gene_biotype = ifelse(!is.na(f$gene_type), f$gene_type, s$gene_type),
  baseMean_fc = f$baseMean, baseMean_salmon = s$baseMean, log2FC_shrunk_fc = f$log2FoldChange_shrunk, log2FC_shrunk_salmon = s$log2FoldChange_shrunk, padj_fc = f$padj, padj_salmon = s$padj)
t$rule1_fc <- ifelse(is.na(f$r1), FALSE, f$r1); t$rule1_salmon <- ifelse(is.na(s$r1), FALSE, s$r1)
t$rule2_pass <- !is.na(t$baseMean_fc) & !is.na(t$baseMean_salmon) & t$baseMean_fc >= 20 & t$baseMean_salmon >= 20
t$rule3a_pass <- t$rule1_fc & t$rule1_salmon & sign(t$log2FC_shrunk_fc) == sign(t$log2FC_shrunk_salmon)
write.table(t, "results_rc1prep_19/candidate_table_rules123a.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
pc <- t$gene_biotype == "protein_coding"
cat("tested_fc", sum(!is.na(fc$padj)), "tested_salmon", sum(!is.na(sa$padj)), "\n")
cat("rule1_fc", sum(t$rule1_fc), "rule1_salmon", sum(t$rule1_salmon), "either", nrow(t), "after_r2", sum(t$rule2_pass), "after_r3a", sum(t$rule2_pass & t$rule3a_pass), "(protein-coding", sum(t$rule2_pass & t$rule3a_pass & pc), ")\n")
