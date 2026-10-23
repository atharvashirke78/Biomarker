# Step 12a v2 summary: lock proposal with candidate features restricted to the 1,390 pool (DRAFT, nothing committed, 15 validation samples not used)

Branch rc1-prep, commit 48ab9fd. Folder results_rc1prep_19/lock_proposal_v2/. All numbers: 19 development samples. Settings as in fit_final_model.R: LASSO lambda.1se, seed 42 (LASSO folds), 43 (cutoff folds), inner 5-fold CV stratified by group x study, Youden J per inner fold, median = final cutoff.

## 1. Final primary model (PROPOSAL)
- **9 genes** (8 up, 1 down) from 1,390 candidate features. lambda.1se = 0.0095792166 (lambda.min would give 8 genes).
- **Final cutoff = 1.352326**; inner-fold cutoffs: 1.35233, 1.2408, 1.35233, 1.35233, 1.35233. The in-sample Youden cutoff is the same value. Score > cutoff = tumour.
- All 9 genes are in the 1,390 pool (by construction).

| Gene | Ensembl ID | Set | LASSO coef | In 1,390 pool | Earlier trial: folds selected (of 19) | Seeds selecting it (of 30) |
|--|--|--|--|--|--|--|
| OSR1 | ENSG00000143867 | down | -1.0008 | yes | 19 / 19 | 30 / 30 |
| KRT80 | ENSG00000167767 | up | 0.6134 | yes | 19 / 19 | 30 / 30 |
| KRT23 | ENSG00000108244 | up | 0.1659 | yes | 18 / 19 | 23 / 30 |
| GRIN2D | ENSG00000105464 | up | 0.4810 | yes | 18 / 19 | 27 / 30 |
| FUT1 | ENSG00000174951 | up | 0.5215 | yes | 18 / 19 | 27 / 30 |
| ETV4 | ENSG00000175832 | up | 0.0094 | yes | 17 / 19 | 26 / 30 |
| LARGE2 | ENSG00000165905 | up | 0.1935 | yes | 17 / 19 | 26 / 30 |
| VSNL1 | ENSG00000163032 | up | 0.0601 | yes | 14 / 19 | 22 / 30 |
| GTF2IRD1 | ENSG00000006704 | up | 0.0545 | yes | 1 / 19 | 22 / 30 |

- Compared with the first draft (candidate features rule 1 + baseMean >= 20, 2,270 genes): MTHFD1L is no longer possible (it is not in the 1,390 pool) and GTF2IRD1 enters with a very small coefficient (0.055). Genes in both versions: ETV4, FUT1, GRIN2D, KRT23, KRT80, LARGE2, OSR1, VSNL1. GTF2IRD1 was selected in only 1 of 19 folds of the earlier trial, so it is a weak member.
- **Development-data performance (OPTIMISTIC):** AUC 1.000; tumour 9/9 called tumour; normal 10/10 called normal. Genes and cutoff were chosen on these samples.

## 2. Seed sensitivity (30 seeds, 1,390 pool)
| Measure | Result |
|--|--|
| Number of genes selected | 2 to 9 (9 genes in 21 of 30 seeds; median 9) |
| Cutoff range | 0.41 to 1.35 (median 1.352) |
| Selection frequency | KRT80 30/30; OSR1 30/30; FUT1 27/30; GRIN2D 27/30; ETV4 26/30; LARGE2 26/30; KRT23 23/30; GTF2IRD1 22/30; VSNL1 22/30 |
| Most common gene sets | 21 seeds: ETV4,FUT1,GRIN2D,GTF2IRD1,KRT23,KRT80,LARGE2,OSR1,VSNL1 ; 3 seeds: ETV4,FUT1,GRIN2D,KRT80,LARGE2,OSR1 ; 3 seeds: KRT80,OSR1|
Table: seed_sensitivity_v2.tsv. Same instability as before: the gene set depends on the seed, the fixed seed 42 gives the most common 9-gene set.

## 3. Cross-study check (information only; same procedure, features = stage-1 fold features of the new fit restricted to the 1,390 pool)
| Direction | n train/test | Features | Up/Down | Cutoff | AUC (95% CI) | Tumour correct | Normal correct |
|--|--|--|--|--|--|--|--|
| A->B | 10/9 | 714 | 8/0 | 3.00 | 1.00 (1.00-1.00) | 4/4 | 3/5 |
| B->A | 9/10 | 567 | 1/2 | -2.34 | 0.96 (0.85-1.00) | 5/5 | 4/5 |

Note: the 1,390 pool was defined using all 19 samples, so this check is not fully independent. The ranking is good (AUC 0.96 to 1.00). A to B: all 4 tumours correct but only 3 of 5 normals; B to A: all 5 tumours and 4 of 5 normals. Against the pass rule (at most 1 error per class) A to B would fail, B to A would pass. This is the cutoff-transfer problem.

## 4. Reference genes within each study (19 development samples; unshrunk DESeq2 log2FC tumour vs normal, ~ group within the study, 90% CI = log2FC +/- 1.645 x SE)
Study SRP219837: 5 normal + 5 tumour; SRP301216: 5 normal + 4 tumour. Pre-filter in each study: >= 10 counts in >= smallest group of that study (5 and 4).
| Gene | SRP219837 log2FC (90% CI) | in +/-0.5 | in +/-0.3 | SRP301216 log2FC (90% CI) | in +/-0.5 | in +/-0.3 | Pooled fit: 0.5 / 0.3 |
|--|--|--|--|--|--|--|--|
| RBM5 | 0.07 (-0.25 to 0.39) | yes | NO | 0.02 (-0.36 to 0.39) | yes | NO | yes / yes |
| TIAL1 | 0.25 (0.10 to 0.40) | yes | NO | -0.15 (-0.51 to 0.22) | NO | NO | yes / yes |
| UBA5 | 0.11 (-0.14 to 0.36) | yes | NO | 0.25 (-0.14 to 0.64) | NO | NO | yes / NO |
| SLC30A9 | 0.11 (-0.10 to 0.31) | yes | NO | -0.19 (-0.65 to 0.28) | NO | NO | yes / yes |
| HNRNPM | 0.20 (-0.11 to 0.51) | NO | NO | 0.30 (-0.01 to 0.61) | NO | NO | yes / NO |
| AKAP8 | -0.10 (-0.43 to 0.23) | yes | NO | 0.17 (-0.23 to 0.56) | NO | NO | yes / yes |
| SNX17 | -0.44 (-0.71 to -0.16) | NO | NO | -0.04 (-0.40 to 0.32) | yes | NO | yes / NO |
| CNBP | -0.35 (-0.75 to 0.05) | NO | NO | -0.03 (-0.38 to 0.31) | yes | NO | yes / NO |
| NMT1 | -0.05 (-0.44 to 0.34) | yes | NO | -0.02 (-0.38 to 0.35) | yes | NO | yes / NO |
| PPP2R5C | -0.08 (-0.42 to 0.27) | yes | NO | 0.22 (-0.26 to 0.71) | NO | NO | yes / yes |

- Inside +/-0.5: SRP219837 7 of 10; SRP301216 4 of 10; both studies: 2 of 10 (RBM5, NMT1).
- Inside +/-0.3: SRP219837 0 of 10; SRP301216 0 of 10. With 4 to 5 samples per group the intervals are wide (width about 0.5 to 0.9 log2), so none can be inside +/-0.3 within a single study.
- baseMean >= 100 on both routes (frozen filter), from the pooled fit:
| Gene | baseMean featureCounts | baseMean Salmon |
|--|--|--|
| RBM5 | 1116 | 1388 |
| TIAL1 | 762 | 776 |
| UBA5 | 354 | 411 |
| SLC30A9 | 718 | 757 |
| HNRNPM | 2021 | 2086 |
| AKAP8 | 468 | 466 |
| SNX17 | 818 | 869 |
| CNBP | 2660 | 2752 |
| NMT1 | 823 | 875 |
| PPP2R5C | 945 | 957 |

## 5. Draft lock files (lock_proposal_v2/)
- final_model_genes.tsv, final_model_spec.txt (genes, signs, formula, cutoff, seeds).
- lock_statement.md: methods, pass rule, QC exclusion rule, conditions, limitations (seed instability, cutoff transfer, wide intervals).
- reference_genes.yaml.draft: new-fit numbers, baseMean >= 100 both-routes filter and margins marked frozen. candidate_rules.yaml.draft: classifier_test filled with methods, pass rule and the QC-exclusion rule; rule 4 stays PLACEHOLDER (no numbers invented). Diffs only: lock_diff_reference_genes.diff, lock_diff_candidate_rules.diff. Tracked config files untouched.
- Code and logs: fit_final_model_v2.R, seed_sensitivity_v2.R, ref_genes_per_study.R and their .log files. glmnet warnings ("class has fewer than 8 observations") are expected.
- sha256 of every file in the folder (except the summary and the checksum file itself): sha256sums.txt. Key files:
  - `final_model_genes.tsv`: 2a510163920b1e7bafef2bf61bfd12203d3c4dbe0da57b18564c297165a5d9d0
  - `final_model_spec.txt`: 02a297517011984ac9ffae74295d6955562eb632d10cbdd014ba33f8658a7f95
  - `reference_gene_list_top10.txt`: 0406c1d3b7d10dd8a251c1d4ab5b117eb33ace45ced0bcea40633c4614f4c6fa
  - `reference_genes.yaml.draft`: 6ebe5f8f34926c29542a9f5a2e0706a0d5bed2c588fb6e9e5fc43122255ac637
  - `candidate_rules.yaml.draft`: 729e568293c4ff02ea2810d263c7836a04024bf3d4e8317f17ba397010304853
