# Nested internal validation, 19 development samples, default settings (2026-10-23)

Run: `scripts/run_nested_default.sh` on pod 2 (32 vCPU, 24 cores used), code commit `54683ba` plus the core-cap change, 395 s wall, exit 0, 525 DESeq2 fits. Files: `docs/results/nested_19dev_20261023/`. Samples: 19 (10 normal, 9 tumour, two studies), featureCounts counts. Settings: 5 outer folds x 5 repeats, 20 label permutations x 5 repeats, seed 1000.

**What this is:** an internal estimate in which DESeq2 and candidate rules 1-2 are repeated inside every training fold, followed by LASSO, score and cutoff. Rules 3a and 3b are not run, so it describes the simplified procedure, not the full locked one (about 1,988 candidates per fold after rules 1-2, against 1,390 after all rules on all samples). It is **not** independent validation.

| Result | Value |
|---|---|
| Pooled out-of-fold AUC, 5 repeats | 1.000 in every repeat |
| Sensitivity (Clopper-Pearson, mean over repeats) | 1.00 (0.63 to 1.00) |
| Specificity | 0.935 (0.58 to 1.00): about one of the ten normals mis-called in some repeats |
| Outer folds run / skipped (observed) | 23 / 2 (LASSO chose no gene at lambda.1se); samples in a skipped fold are not scored |
| Mean genes selected per fold | 7.1 |
| Permutation null, 20 of 20 usable | null AUC mean 0.40, 95th percentile 0.73; empirical p = 0.048 |

**Gene stability** (`nested_selection_frequency.tsv`; 23 folds): 43 genes were selected at least once. All 9 locked-model genes appear, but only four in at least half of the folds: OSR1 100%, KRT80 96%, GRIN2D 74%, FUT1 61%. LARGE2 39%, ETV4 35%, KRT23 26%, VSNL1 17%, GTF2IRD1 4% (1 of 23 folds). MTHFD1L (35%) and CRIPTO (30%) are selected as often as some locked genes.

## How to read it
1. The separation between tumour and normal in this cohort is strong and survives repeating the gene screen inside the folds: AUC 1.0 in all repeats. That is a better statement than the in-sample AUC of 1.0 of the development fit, but the cohort is small (19 samples, two studies), the screen is simplified, and a perfect AUC with wide interval bounds still needs an independent cohort.
2. **The permutation p-value 0.048 is the smallest this design can give** (B = 20 gives 1/21): it says the observed AUC beat all 20 permuted runs, nothing finer. The null rests on few folds: 324 of the 500 permuted folds (176 ran) had fewer than 2 candidates after rules 1-2 or no LASSO gene and were skipped, so each null AUC rests on 2 to 5 repeats and often on only a few folds. A null mean below 0.5 (0.40) is commonly seen with cross-validation on shuffled labels in small samples and is not by itself evidence of leakage (not verified here against a source). To tighten it, run B of 100 to 200 (about 30 to 60 minutes on pod 2, an estimate).
3. The locked 9-gene panel is only partly stable under resampling (a core of four). The panel in `lock/diagnostic_v1/` is unchanged; this is information for the validation plan and the report text, not a reason to change the lock.
4. Both tumour and normal samples come from the same two studies as the genes were screened on; a leave-one-study-out evaluation was not done (not possible with two studies and the group-by-study design).

## Not done or not verified
- Rules 3a and 3b inside the folds; patient-aware folds (no repeated patient_id here); one-batch design; identical output at different core counts (the earlier small test used 2 cores, this run 24: the full-size permutation null was not re-run to compare).
- Independent validation on the 15 held-out samples.
