# Threshold inventory (read-only snapshot)

Values read from the files at the lines given (nextflow/ paths: `modules/`, `conf/`, `nextflow/bin/`). Nothing was changed. Candidate-rule entries are DRAFT. [Updated 2026-10-15: rules 1-3b are DECIDED and rule 4, the classifier and the reference genes are LOCKED at tag lock-diagnostic-v1; rows below re-read from config/candidate_rules.yaml and config/reference_genes.yaml.]

## 1. Thresholds the pipeline applies

| Step | Parameter | Value | File:line |
|---|---|---|---|
| fastp | `-q` base qualified if Phred >= | 20 | `conf/thresholds.config:6`, used `modules/fastp.nf:24` |
| fastp | `--unqualified_percent_limit` | 40 | `thresholds.config:7`, `fastp.nf:24` |
| fastp | `--cut_tail`, window / mean quality | 4 / 20 | `thresholds.config:8-9`, `fastp.nf:25` |
| fastp | `--length_required` | 36 | `thresholds.config:10`, `fastp.nf:26` |
| fastp | `--n_base_limit` | 5 | `thresholds.config:11`, `fastp.nf:26` |
| fastp | adapter detection | `--detect_adapter_for_pe` (no number) | `fastp.nf:23` |
| STAR | `--outFilterMultimapNmax` | 20 | `thresholds.config:14`, `modules/star.nf:55` |
| STAR | all other alignment options | STAR defaults (none set) | `star.nf:53-57` |
| featureCounts | `-Q` min MAPQ | 10 | `thresholds.config:15`, `modules/featurecounts.nf:27` |
| featureCounts | other flags | `-p --countReadPairs -B -C`, `-t exon -g gene_id`, `-s` from Salmon library type | `featurecounts.nf:19-28` |
| Salmon | options | `-l A --validateMappings --gcBias --seqBias`, no numeric cutoff | `modules/salmon.nf:18-19` |
| Memory guard | `min_mem_star_gb` | 32 | `thresholds.config:16` |
| Per-sample QC (WARN) | STAR unique % below | 70 | `thresholds.config:19`, `bin/sample_qc.py:36` |
| Per-sample QC (WARN) | featureCounts assigned % below | 60 | `thresholds.config:20`, `sample_qc.py:36` |
| Per-sample QC (WARN) | Salmon mapped % below | 70 | `thresholds.config:21`, `sample_qc.py:36` |
| Per-sample QC (FAIL, excluded) | STAR unique % below | 50 | `thresholds.config:22`, `sample_qc.py:37` |
| Per-sample QC (FAIL, excluded) | Salmon mapped % below | 50 | `thresholds.config:23`, `sample_qc.py:37` |
| Per-sample QC | fastp Q30 and passed-read % | calculated and reported, no cutoff | `sample_qc.py:24-28` |
| Cohort QC pre-filter | >= 10 counts in >= min(3, n) samples, for VST/PCA | 10 / 3 | `thresholds.config:26-27`, `bin/cohort_qc.R:14` |
| Cohort QC outlier flag | abs(robust z, MAD) on PC1 or PC2 within each study x group cell | 3 (flag only, nothing removed) | `thresholds.config:28`, `cohort_qc.R:30` |
| ANALYSIS input checks (PROVISIONAL) | samples per group / % unmapped gene IDs | >= 3 / <= 5 % | `thresholds.config:31-32`, `validate_analysis_input.py:85,118` |
| DESeq2 pre-filter | >= 10 counts in >= n samples, n = smallest group in the sample sheet (`de_min_samples = 'auto'`; 9 for the 19 samples, VERIFIED: the 2026-10-14 run printed "pre-filter: >= 10 counts in >= 9 samples", `results_rc1prep_19/run_fc.log`) | 10 / auto | `thresholds.config:35-36`, `bin/deseq2.R:21` (rc1-prep, changed from 5) |
| DESeq2 test | `results(alpha = alpha)`, same alpha as the padj cut (`de_padj`) | 0.05 | `deseq2.R:8,27` (rc1-prep, was package default 0.1) |
| DESeq2 significance | padj / abs(log2FC) | 0.05 / 1, unshrunk value (Step 4 to switch to shrunk) | `thresholds.config:37-38`, `deseq2.R:33` |
| Shrinkage | `lfcShrink(type = "apeglm")` | no numeric cutoff | `deseq2.R:27` |
| ORA | significant set padj / abs(LFC); universe = tested genes with a symbol | 0.05 / 1 | `bin/ora.R:7-10` |
| ORA | gene-set size min / max; BH q and p cutoff | 10 / 500; 0.05 | `thresholds.config:43-45`, `ora.R:19,22` |
| ORA | KEGG source | live KEGG web service, fallback MSigDB KEGG sets | `ora.R:14,22-26` |
| GSEA | ranking: Wald `stat`, one row per symbol (largest abs(stat)) | - | `bin/gsea.R:7-8` |
| GSEA | `fgseaMultilevel`, `nPermSimple`, seed, padj | 10,000, 42, 0.05 | `thresholds.config:46-47`, `gsea.R:11,14,17` |
| GSEA | gene-set size min / max | 10 / 500 | `thresholds.config:43-44`, `gsea.R:14` |
| GSEA | collections | Hallmark, KEGG, GO BP | `gsea.R:10` |
| Candidate rules | rule 1 padj < 0.05 and abs(shrunk LFC) >= 1 (DECIDED): `padj_max: 0.05`, `abs_shrunk_log2fc_min: 1.0` | 0.05 / 1.0 | `config/candidate_rules.yaml:12-13` |
| Candidate rules | rule 2 baseMean >= 20 on both routes (DECIDED): `basemean_min: 20`; baseMean >= 50 on both is a reported tier only, NOT a filter (`stricter_tier: basemean_min: 50`) | 20 (filter) / 50 (flag) | `candidate_rules.yaml:17`, `:19-20` |
| Candidate rules | rule 3b: rule 1 holds in >= `loo_min_pass` = 18 of 19 leave-one-out refits (DECIDED 2026-10-14; 17/19 = 0.895 is below 90%, 18/19 = 0.947) | 18 | `candidate_rules.yaml:32` |
| Candidate rules | rule 4 (independent reproduction): `status: LOCKED`, same direction and padj < 0.05 in the 15 samples, `pass_threshold: none (reported only)` | - | `candidate_rules.yaml:34-39` |
| Classifier test | `status: LOCKED`; balanced up/down score, LASSO lambda.1se, seeds 42 (LASSO folds) and 43 (cutoff folds), inner 5-fold CV, Youden median cutoff; pass rule at most 1 error per class | cutoff 1.352325902 (`lock/diagnostic_v1/final_model_spec.txt`) | `candidate_rules.yaml:40-48` |
| Reference genes | `basemean_min_both_routes: 100` (FROZEN); `equivalence_margins_log2fc: [0.5, 0.3]` (FROZEN); 90 % CI | 100; 0.5 / 0.3 | `config/reference_genes.yaml:10-11`, `:6`; applied (hard-coded copy) `scripts/downstream/reference_gene_check.R:13` |

Notes:
- The DESeq2 significance rule in the shipped script still uses the unshrunk LFC; the shrunk rule applies only to the candidate rules so far.
- FIXED on branch rc1-prep (commit ecedc96) and run on 2026-10-14 (`results_rc1prep_19/deseq2_summary.json`: padj_cut 0.05): `results()` now gets the same alpha as the padj cut. Line numbers above are for rc1-prep (checked 2026-10-15).
- Per-sample QC has no minimum read count and no minimum trimmed-read count.

## 2. Check against the linked pages

| Page | Result |
|---|---|
| fastp (GitHub README) | Opened. `-q` default 15 ("phred quality >=Q15 is qualified"); `-u` default 40; `-l` default 15; `-n` default 5. Quality filtering, length filtering and adapter trimming on by default. |
| nf-core/rnaseq parameters | Opened. `min_trimmed_reads` default 10000 ("Minimum number of trimmed reads below which samples are removed from further processing"). `min_mapped_reads` default 5 ("Minimum percentage of uniquely mapped reads below which samples are removed from further processing"). `stranded_threshold` 0.8, `unstranded_threshold` 0.1. `extra_star_align_args` and `extra_salmon_quant_args` have no defaults on the page. |
| DESeq2 vignette | Opened, only the first 100,000 of 432,967 characters read. Pre-filtering: "keep only rows that have a count of at least 10 for a minimal number of samples"; example `smallestGroupSize <- 3; keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize`; set the minimal number to the smallest group size. `results()` does independent filtering automatically, `alpha` default 0.1. apeglm is the shrinkage method named. |
| ENCODE long RNA-seq STAR options | NOT FOUND. My first URL was a guessed document link and returned 404. `/rna-seq/long-rnas/` and `/data-standards/rna-seq/long-rnas/` opened but list no STAR options and no mapping-rate standard. They state only read depth: 30 million aligned reads per replicate (standard bulk RNA-seq), 10 million (knockdown experiments). STAR parameter list unverified. |

Comparison with our values:

| Item | Ours | Page |
|---|---|---|
| fastp `-q` | 20 | 15 default |
| fastp `--length_required` | 36 | 15 default |
| fastp `-u` and `-n` | 40 and 5 | equal to the defaults |
| STAR unique FAIL | 50 % | nf-core `min_mapped_reads` 5 % |
| Minimum trimmed reads | none | nf-core 10000 |
| DESeq2 pre-filter | >= 10 counts in >= 9 samples (rc1-prep; was >= 5 before 2026-10-14) | vignette example: smallest group size; our smallest group (tumour) has 9 samples |
