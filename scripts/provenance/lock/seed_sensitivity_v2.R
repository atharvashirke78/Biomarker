# DRAFT v2 info only (features = 1,390 pool): how much does the final gene set / cutoff depend on the seed? (19 development samples only; same procedure as fit_final_model.R)
src <- readLines("results_rc1prep_19/lock_proposal_v2/fit_final_model_v2.R"); eval(parse(text = src[1:grep("^# ---- final model on all 19", src) - 1]))
all19 <- colnames(cn); strata <- paste(y[all19], study[all19]); sym <- setNames(res$gene_name, res$gene_id)
out_rows <- list(); sets <- list()
for (sd in 1:30) { fm <- suppressWarnings(fit_model(feat19, all19, strata, seed = sd)); sets[[sd]] <- sym[names(fm$cf)]
  out_rows[[sd]] <- data.frame(seed = sd, n_genes = length(fm$cf), n_up = length(fm$up), n_down = length(fm$dn), cutoff = fm$thr, genes = paste(sort(sym[names(fm$cf)]), collapse = ",")) }
S <- do.call(rbind, out_rows); write.table(S, "results_rc1prep_19/lock_proposal_v2/seed_sensitivity_v2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
tb <- sort(table(unlist(sets)), decreasing = TRUE); print(tb); print(summary(S$n_genes)); print(summary(S$cutoff)); cat("'RUVBL1' in candidate features:", "ENSG00000175792" %in% sub("\\..*", "", feat19), "\n")
