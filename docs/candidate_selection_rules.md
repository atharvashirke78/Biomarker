# Candidate biomarker selection rules (APPROVED for rules 1-3)

Scope: 19 diagnostic samples (9 tumour, 10 normal; studies A 5+5, B 4+5), DESeq2 `~ study + group`, reference normal. Nothing here touches the 15 validation samples. Companion file: `config/candidate_rules.yaml`. Tables are keyed on versioned `gene_id`; the symbol is a label only (PDE4C = ENSG00000105650.23 and ENSG00000285188.1). A `gene_biotype` column is kept for all genes; protein-coding is NOT restricted (decision at the final-panel stage).

Status words: **DECIDED** = owner decision; **PLACEHOLDER** = not run now. (Earlier draft with PROPOSED values: `*.bak_draft`.)

## Rule 1: significance (DECIDED)
padj < 0.05 and |apeglm-shrunk log2FC| >= 1. The candidate table holds every gene that passes rule 1 on featureCounts or on Salmon (one row per gene_id).

## Rule 2: measurable (DECIDED): baseMean >= 20 on both routes
Cutoff table (rule-1 genes removed) from Step B: at 20, 304 (featureCounts) and 210 (Salmon). Stricter tier: column `baseMean_ge_50_both` (TRUE/FALSE, baseMean >= 50 on both routes) is reported for the final panel and is NOT a filter.

## Rule 3: robust
(a) DECIDED: rule 1 passes on both featureCounts and Salmon, same direction.
(b) DECIDED: leave-one-out on featureCounts only. Refit DESeq2 + apeglm 19 times, each time dropping one of the 19 samples, with the same design as the 19-sample run (`~ batch + group`, batch = study for this cohort), the same pre-filter (>= 10 counts in >= n samples, n = size of the smallest group of the FULL 19-sample cohort = 9, calculated once and used in all 19 refits) and the same `results(alpha = 0.05)` / `lfcShrink(type = "apeglm")` call as the 19-sample run. alpha (`de_padj`), the count 10 (`de_min_count`) and n (`de_min_samples`) are read from `conf/thresholds.config`, the same parameters `deseq2.R` uses, so they cannot drift apart. (Changed on branch rc1-prep; the earlier run used pre-filter 5 and the default alpha 0.1.) A gene passes a refit if rule 1 holds in the same direction as the 19-sample run (a gene not tested in a refit fails it). Stable = passes in at least `loo_min_pass` of 19 refits (parameter in `config/candidate_rules.yaml`; **18** of 19, DECIDED 2026-10-14). Per gene we record the fraction and which sample removals broke it.
Note: the owner's text said "at least 90% (17 of 19 or more)". 17/19 = 0.895 is below 90%, so on 2026-10-14 the owner set the threshold to 18 of 19 (0.947), which meets "at least 90%". Earlier runs used 17.

## Rule 4 (PLACEHOLDER): reproduce in the 15 independent samples
Pass numbers to be fixed and committed before the 15 samples are touched.

## Classifier test (PLACEHOLDER): AUC, sensitivity, specificity
Pass numbers to be fixed and committed before the 15 samples are touched.
