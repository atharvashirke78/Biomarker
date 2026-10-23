# Errata, diagnostic v1 (tag lock-diagnostic-v1)

Written 2026-10-15, after the lock; item (a) added and items relettered on 2026-10-15. These are corrections to records, wording and labels only. The locked gene list, signs, formula, cutoff (1.352325902), seeds, rules and pass rule do not change, and no sample call changes. The tag is not moved; the locked files are not edited.

**(a) No recorded command for the 2026-10-14 DESeq2 re-run.** `docs/run_log.md` has no entry giving the command that produced the new 19-sample DESeq2 fit in `results_rc1prep_19/`: `deseq2_results_all.tsv` at 2026-10-14 05:29:59 and `deseq2_salmon_results_all.tsv` at 05:30:30, with logs `run_fc.log` and `run_sal.log`. The exact command was searched for and NOT FOUND in any of these:
- shell history (no `.bash_history` on this pod);
- `logs/`, `scripts/`, `docs/` and the `results_rc1prep_19/*.md` notes;
- the Claude session transcripts in `/workspace/.claude` (the only earlier transcript ends on 2026-10-08).

What the files do show:
- `deseq2.R` was committed with the alpha and pre-filter fix at 05:21:14 (commit ecedc96).
- `sample_sheet_rc1prep.tsv` was written at 05:29:33.
- `run_fc.log` reports "pre-filter: >= 10 counts in >= 9 samples: kept 15111 of 78670 genes".
- `deseq2_summary.json` records padj_cut 0.05 and lfc_cut 1.

A probable form, RECONSTRUCTED from the usage line of `nextflow/bin/deseq2.R` and the output names, NOT VERIFIED as the command actually run:
```
cd results_rc1prep_19
../env/prefix/r_dge/bin/Rscript ../nextflow/bin/deseq2.R matrix   ../results/04_counts/counts_featurecounts.tsv sample_sheet_rc1prep.tsv ../ref/gene_id_to_symbol.tsv deseq2        10 auto 0.05 1 15 50 > run_fc.log  2>&1
../env/prefix/r_dge/bin/Rscript ../nextflow/bin/deseq2.R tximport ../results/04_counts/salmon_tximport.rds    sample_sheet_rc1prep.tsv ../ref/gene_id_to_symbol.tsv deseq2_salmon 10 auto 0.05 1 15 50 > run_sal.log 2>&1
```
The result files are consistent with the locked rules: the 19-sample refit inside the rule 3b run reproduces them (max abs difference in log2FC 4.1e-14; `results_rc1prep_19/rule3b/reproduction_check_19_sample_fit.tsv`).

**(b) Lowest tumour score in `lock_statement.md`, limitation 3.** The statement says "lowest tumour score 3.67". The locked model's score file (`results_rc1prep_19/lock_proposal_v2/development_scores.tsv`) gives 3.541, for SRR10043502. The value 3.67 (3.672) belongs to the earlier proposal (`results_rc1prep_19/lock_proposal/development_scores.tsv`). The highest normal score (-0.84) is correct. The scores were recomputed from `results/04_counts/counts_featurecounts.tsv` with the locked genes, signs and formula, without refitting. They match the score file (largest difference 4.9e-15), and all 19 calls are the same. No cutoff or call changes.

**(c) First line of `final_model_spec.txt`.** Line 1 reads "# DRAFT final primary model (not locked)". That wording is left over from the proposal. The file is the locked version at tag lock-diagnostic-v1 (sha256 2a5101… for the genes file and 02a297… for this file, as in `sha256sums_final.txt`). LOCK_NOTES.md (b) covers the draft wording in `lock_statement.md` and `step12a_v2_summary.md` but did not name this line.

**(d) Rule 3b funnel label.** In `results_rc1prep_19/rule3b_rerun18/funnel_counts.tsv` and `results_rc1prep_19/rule3b_rerun18_cfg/funnel_counts.tsv`, the row label reads "+ rule 3b (LOO >= 17 of 19, featureCounts)". Its count of 1,390 (1,353 protein-coding) is the 18-of-19 result, as locked (`config/candidate_rules.yaml` loo_min_pass: 18). The 17-of-19 result is 1,462 (`results_rc1prep_19/rule3b/funnel_counts.tsv`). The label is wrong; the count is right.

**(e) Two study A runs left out without a documented project reason.** ENA lists 12 runs for SRP219837 (7 tumour, 5 normal). Two tumour runs are not in the development set: SRR10043503 (COAD-T015) and SRR10043504 (COAD-T016). They are not in `samplesheets/diagnostic_20.csv`. They are also not in the paper's run list: Table S2 of `docs/reference/hosseini_2023_MOESM2_ESM.xlsx` contains our 20 runs and neither of these two. That is why they were left out. No other project reason (QC, data quality or sample type) is documented. Including them would have given study A 7 tumours and 5 normals. Nothing in the lock depends on this note.

**(f) Audit findings, 2026-10-15.** Full audit: `results_rc1prep_19/audit_compare_20261015/audit_report.md` (git-ignored, on the pod only). These findings are about records, documentation and code reuse. None of them changes a locked value.
1. **Six inputs of the diagnostic report have no logged command in `docs/run_log.md`:**
   - the 2026-10-14 DESeq2 fit (`results_rc1prep_19/deseq2_results_all.tsv`, `deseq2_salmon_results_all.tsv`, their summaries and the 6 DE plots; see item (a));
   - `results_rc1prep_19/sample_sheet_rc1prep.tsv`;
   - the step that created `results_rc1prep_19/final_candidates/final_candidates_18of19.tsv` (only its identity with `rule3b_rerun18/` is logged);
   - `results_rc1prep_19/lock_proposal_v2/` (`fit_final_model_v2.R`, `development_scores.tsv`, `seed_sensitivity_v2`, `ref_genes_per_study`);
   - the 2026-10-09 09:18 re-run of the cohort outlier flags (`results/05_cohort_qc/cohort_qc_flags.tsv`, `cohort_qc_summary.json`);
   - `results_rc1prep_19/biomarker_19_samples_all_data.xlsx`.
2. **`docs/threshold_inventory.md` was stale in 6 places at the time of the audit:**
   - rule 2 listed as "baseMean >= 50 (PROPOSED)" (it is 20; 50 is a flag tier);
   - rule 1 line reference;
   - rule 4 and the classifier listed as "placeholders" (both LOCKED);
   - pre-filter n = 9 marked "UNVERIFIED";
   - the comparison table pre-filter of ">= 5 samples";
   - `deseq2.R` line references and the "not yet run" note.
3. **Tool-version docs:** glmnet 4.1.10 was missing from `docs/tool_versions.tsv`, and `docs/run_log.md:23` lists the pre-build pins (fastp 0.23.4, DESeq2 1.42.0, apeglm 1.24.0, R 4.3, subread 2.0.6), which were superseded by the built versions in `run_log.md:93`.
4. **Thresholds not read from config:**
   - `scripts/downstream/reference_gene_check.R:13` hard-codes baseMean 100, margins 0.5 / 0.3 and CI 0.90 (equal to `config/reference_genes.yaml:10-11`; comments at :5 and :42 still say PROPOSED);
   - `scripts/downstream/reactome_analysis.R:14` hard-codes set size 10-500, q 0.05, 10,000 permutations and seed 42 (equal to `conf/thresholds.config:43-47`);
   - `scripts/candidate_selection.R:104` hard-codes the funnel label "LOO >= 17 of 19", the source of item (d); the threshold used, 18, is read from config.
5. **Downstream scripts are tied to the 19 development samples.** `scripts/candidate_selection.R`, `scripts/downstream/lasso_stage1_deseq2.R`, `lasso_stage2_glmnet.R`, `lasso_arm3.R`, `lasso_arm4.R`, `reference_gene_check.R` and `lock_proposal_v2/fit_final_model_v2.R` hard-code:
   - 19 samples (`ncol == 19`);
   - the two study IDs (SRP219837, SRP301216);
   - the group names "tumour" / "normal" (coefficient `group_tumour_vs_normal`).

   They would stop or give wrong results on another sample sheet.
