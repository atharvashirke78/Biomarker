# Fully nested internal validation: design, limits, cost

Script: `nextflow/bin/nested_full_validation.R`. Test: `tests/regression/test_nested_full_validation.sh`.
Status: additive and report-only. It reads the same inputs as `internal_validation.R`, writes only its own four files, never changes the model and never fails a pipeline.

## Why

`internal_validation.R` repeats LASSO + score + cutoff inside the resamples but keeps the candidate screen (DESeq2 on all samples, rules 1-3b) outside, so its numbers are optimistic (its own LIMITATION text says so). Selection and tuning must be repeated inside every training fold (Ambroise & McLachlan 2002; Varma & Simon 2006; MAQC-II; nestedcv; reviewed in `docs/review/02_biomarker_methods.md`). A label-permutation null shows what the whole procedure produces when there is no signal.

## Procedure (per repeat and outer fold, training samples only)

1. Folds: outer K-fold (default 5) stratified by group x batch; if `patient_id` repeats, all samples of a patient stay together (same helper style as `fit_diagnostic_model.R`, `--patient-level-cv true` default here; `false` gives sample-level folds).
2. DESeq2 on the training samples: design `~ batch + group` (`~ group` with one batch), pre-filter `de_min_count` in `de_min_samples` samples, `results(alpha = de_padj)`, `lfcShrink(type = "apeglm")` (values from `conf/thresholds.config`, read, never edited).
3. Rule 1: `padj < padj_max` and `|shrunk LFC| >= abs_shrunk_log2fc_min`; rule 2: `baseMean >= basemean_min` (single route), values from `config/candidate_rules.yaml`, same expressions as `candidate_rules.R`.
4. LASSO: glmnet binomial, alpha 1, `foldid` stratified by group x batch (seed `lasso_seed`, K `cv_folds` from `config/questions.yaml`), deviance, `lambda.1se`; log2(CPM+1); up/down score = mean(up) - mean(down); cutoff = median of K inner-fold Youden thresholds (seed `cutoff_seed`).
5. Predict the held-out fold with the training-fold genes and cutoff.
6. Skip the fold (logged reason, never a crash): empty batch x group cell in training, a class with < 2 training samples, < 2 candidates after rules 1-2, no gene at `lambda.1se`, or any fit error.

## Not run (stated in every output)

Rule 3a needs the second quantification route; rule 3b is leave-one-out inside the fold (N more DESeq2 fits per fold). Neither is run. The estimate is for the SIMPLIFIED procedure (rules 1-2 + LASSO + score + cutoff), not the locked procedure; the locked procedure has a stricter candidate list, so its true out-of-sample behaviour could be better or worse. The 19-sample run gives about 2000 candidates per fold, versus 1390 after rules 1-3b on all samples: the screens are not equivalent.

## Outputs (deterministic: no timestamps, fixed seeds)

- `nested_validation_summary.tsv` (metric, value, ci_low, ci_high, label; every label ends with `NESTED, rules 1-2 only; NOT independent validation`): nested AUC (mean over repeats of the pooled out-of-fold AUC; interval = mean of per-repeat DeLong bounds), repeat spread (min/max, not a CI), sensitivity and specificity (pooled calls at each fold's own cutoff; Clopper-Pearson, mean of per-repeat bounds), folds run/skipped, mean selected genes, mean candidates, per-repeat AUCs, null AUC mean and 95th percentile, empirical p.
- `nested_selection_frequency.tsv`: per gene, folds in which the LASSO selected it, frequency, frequency in the candidate list, majority direction (observed run only).
- `nested_permutation_null.tsv`: one row per permutation.
- `NESTED_VALIDATION.md`: plain-language summary, skip reasons, LIMITATION.
- Any analysis problem gives a `NOT RUN: reason` row and exit 0. Usage errors (missing/unknown option, non-integer value) exit 2. `--check-only true` validates inputs and both R environments (exit 1 if unusable).

## Permutation null

Labels are permuted within batch (patient level when a patient_id repeats and every patient has one group and batch), seeded per permutation; the whole procedure (folds, DESeq2, rules, LASSO, cutoff, pooling) is repeated with `--perm-repeats` repeats (default = `--repeats`). Null statistic = mean over its repeats of the pooled AUC, the same statistic as the observed value. `p = (b+1)/(B+1)`, B = permutations with at least one scored fold, so the smallest possible p is 1/(B+1) (0.048 for B = 20). A permutation in which every fold is skipped (no candidates) has no AUC and is excluded from B; this removes "no discrimination" runs from the null and therefore makes p larger (conservative); the count of usable permutations is stated in the output.

## Environments

No environment has all four packages: `r_dge` has DESeq2 and apeglm but not glmnet/pROC; `r_classifier` the reverse (checked with `requireNamespace`). The script is therefore an orchestrator (base R only, can be started by any Rscript) that calls itself twice: `--stage screen` under `r_dge` (DESeq2 + rules 1-2, one RDS per fold in `<out>/nested_work_<hash>/`), then `--stage model` under `r_classifier`. The folds and permuted labels are rebuilt identically in both stages from base-R seeded code. The work folder is a resume cache (keyed on input hashes and settings) and is removed after a successful run unless `--keep-work true`. Override interpreters with `--rscript-dge`, `--rscript-clf`.

## Cost (shared pod; load average about 80 at test time, `--cores 2`)

Dominated by DESeq2 + apeglm: one fit per outer fold per repeat per (permutation + 1). Measured: 20 fits (1 repeat, 3 permutations, 19 samples) took about 138 s wall in total, about 7 s wall per fit at 2 cores. Default settings are (1 + 20) x 5 x 5 = 525 fits, extrapolated to about 60 minutes at the same load (NOT measured). The model stage takes seconds. Cost grows with samples and genes; `--cores` is capped at 4. Reduce with `--perm-repeats 1` (105 fits, about 12 min) if needed; this changes the null statistic's precision, not its meaning.

## Limits

- Few samples: wide, unstable numbers; a fold's test set is 3-5 samples.
- Pooled out-of-fold scores come from different models and cutoffs; DeLong CI treats them as one model (approximate) and repeats reuse the same samples (repeat spread and Clopper-Pearson means are descriptive).
- Rules 3a/3b not run (see above). Batch effects are handled only by the design term, as in the pipeline.
- With a strict rule 1 and a small training fold, few candidates may pass; folds are then skipped, which also happens under the null.
- `--screen-outside true` is a COMPARATOR (screen once on all samples, then nested LASSO): optimistic by design, labelled as such, and only used to demonstrate selection bias.
- Not independent validation.
