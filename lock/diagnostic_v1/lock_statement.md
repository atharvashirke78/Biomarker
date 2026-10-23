# Lock statement, diagnostic v1 (locked by git tag lock-diagnostic-v1)
Branch rc1-prep, base commit 48ab9fd. Nothing here has used the 15 validation samples. The lock is recorded in config/candidate_rules.yaml and config/reference_genes.yaml, committed together with this file.

## Methods
- **Primary method:** balanced up/down score (Arm 4 design). Candidate genes are the 1,390 final candidates (rules 1, 2, 3a, 3b at 18 of 19). Genes are chosen by LASSO (glmnet, binomial, alpha 1, lambda.1se = 0.0095792166, seed 42 for the CV folds) on the 19 development samples only. score = mean(log2(CPM+1) of up genes) - mean(log2(CPM+1) of down genes). A sample is called tumour if score > cutoff. Cutoff = 1.35233: inner 5-fold CV stratified by group and study, Youden J per inner fold, median of the five fold cutoffs, seed 43. Genes, signs and cutoff: final_model_genes.tsv and final_model_spec.txt.
- **Secondary method:** reference-gene score (Arm 3 design, the 10 reference genes below), reported for comparison only; not part of the pass test. A final Arm 3 model was not fitted in this draft.
- **Reference genes (frozen):** RBM5, TIAL1, UBA5, SLC30A9, HNRNPM, AKAP8, SNX17, CNBP, NMT1, PPP2R5C. Filter: baseMean >= 100 on both routes (frozen). Equivalence margins 0.5 and 0.3 log2FC on the 90 % CI of tumour vs normal (frozen). The list, filter and margins are not revised after the test.

## Pass rule for the 15 samples (to be applied once; not run now)
At most 1 classification error per class (tumour, normal), and the same within each study if both studies are present. AUC is reported with its confidence interval but is not the pass test. Confidence intervals will be wide because the test set is small.

## QC exclusion rule
A sample that fails the existing QC flags (sample QC status FAIL: STAR uniquely mapped < 50 % or Salmon mapped < 50 %) is excluded automatically and reported with its reason. WARN samples are kept and reported. The cohort outlier flag (MAD) is flag-only. No sample is excluded or kept by manual choice.

## Conditions
- The 15 validation samples are not used for selection or tuning of any gene, cutoff, list, filter or threshold.
- No change to the gene list, signs, cutoff, reference-gene list, filter, margins, QC rule or pass rule after the test.
- If the pass rule fails, the failure is reported as it is and nothing is adjusted.
- Rule 4 (candidate reproduction in the 15 samples): a candidate gene is "reproduced" if, in the 15 independent samples, it has the same direction as in the 19 development samples (sign of the shrunk log2FC) and padj < 0.05 in a DESeq2 run with the same model as the 19-sample run (~ study + group if more than one study is present, ~ group otherwise; same pre-filter rule by smallest group size). The result is reported as the number and percentage of the 1,390 candidates reproduced, with no pass or fail threshold. Status LOCKED in config/candidate_rules.yaml.

## Known limitations (stated in plain words)
1. **The selected gene set depends on the random seed.** With the same 19 samples and other seeds, the LASSO picks between 2 and 9 genes; the same 9 genes came out in 21 of 30 seeds. The lock fixes the seed (42), so the model is reproducible, but the gene list is one of several that fit the development data equally well. The core genes (KRT80, OSR1) were selected in all 30 seeds.
2. **The cutoff does not transfer well between studies.** When the model is trained on one study and tested on the other, the ranking is good (AUC 0.96 to 1.00) but the cutoff learnt in one study calls some normal samples of the other study tumour (A to B: 3 of 5 normals correct; B to A: 4 of 5). The pass rule allows at most 1 error per class, so this problem could cause a failure on the 15 samples. The earlier trial showed the same pattern.
3. **Development performance is optimistic.** Genes and cutoff were chosen on the same 19 samples. The two classes are far apart there (lowest tumour score 3.67, highest normal score -0.84), so the cutoff is weakly determined by the data.
4. **Confidence intervals will be wide** because the test set is small (15 samples). AUC is not the pass test.
5. **Reference genes within one study:** with 4 to 5 samples per group, the 90 % intervals are wide; no gene is inside +/-0.3 within either study (table in the summary). The freeze rests on the pooled 19-sample fit.

## sha256 of the draft files at preparation time (the committed config files and this statement differ from these drafts; the git tag holds the final versions)
- `results_rc1prep_19/lock_proposal_v2/final_model_genes.tsv`: 2a510163920b1e7bafef2bf61bfd12203d3c4dbe0da57b18564c297165a5d9d0
- `results_rc1prep_19/lock_proposal_v2/final_model_spec.txt`: 02a297517011984ac9ffae74295d6955562eb632d10cbdd014ba33f8658a7f95
- `results_rc1prep_19/lock_proposal_v2/reference_gene_list_top10.txt`: 0406c1d3b7d10dd8a251c1d4ab5b117eb33ace45ced0bcea40633c4614f4c6fa
- `results_rc1prep_19/lock_proposal_v2/reference_genes.yaml.draft`: 6ebe5f8f34926c29542a9f5a2e0706a0d5bed2c588fb6e9e5fc43122255ac637
- `results_rc1prep_19/lock_proposal_v2/candidate_rules.yaml.draft`: 729e568293c4ff02ea2810d263c7836a04024bf3d4e8317f17ba397010304853

Locked 2026-10-14. Code commit at lock time: 48ab9fd. The 15 validation samples were not used to select genes, tune the cutoff, or modify the model.
