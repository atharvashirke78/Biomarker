# Pipeline specification: diagnostic RNA-seq biomarker pipeline (FASTQ to report)

Written 2026-10-16 (documentation only; nothing was run to write it). Repository `/workspace/Biomarker_pipeline`, branch `rc1-prep`, HEAD `25c95ac`, `VERSION` 0.9.0-dev. The locked model is tag `lock-diagnostic-v1` (commit `10bd636`), files in `lock/diagnostic_v1/`.

**Scope and status in one paragraph.** This is research software for tumour-vs-normal RNA-seq (colorectal data so far). Everything is **technical / development**: none of it is clinical validation. Only the pre-refactor version of QUANT has run on real FASTQ (20 samples, 2026-10-08/09). The refactored QUANT has only been run with `-stub-run`. The 15 held-out samples have not been quantified, so the locked 9-gene model has **no independent validation yet**.

**How statuses are labelled.** Every "Check" is labelled VERIFIED (a file or log in this repository shows it, named in the cell) or NOT VERIFIED. "Real data" says whether the step has run on real data in its current (refactored) code: **YES**, **NO (stub only)**, **NO**, or **MANUAL** (a script run by hand, not by Nextflow).

**Where things live**

| What | Where |
|--|--|
| Code | `nextflow/main.nf`, `modules/*.nf`, `nextflow/bin/*`; manual scripts in `scripts/`, `scripts/downstream/` |
| All thresholds | `conf/thresholds.config` (included from `nextflow/nextflow.config:49`); rules in `config/candidate_rules.yaml`; DISCOVER settings in `config/questions.yaml` |
| **The lock** (never edited, never written by any entry) | `lock/diagnostic_v1/` (`final_model_genes.tsv`, `final_model_spec.txt`, `lock_statement.md`, `sha256sums_final.txt`, ...), plus `config/candidate_rules.yaml` and `config/reference_genes.yaml`; hashes in `lock/diagnostic_v1/sha256sums_final.txt` (7 of 7 OK on 2026-10-16: VERIFIED) |
| Real quantification of the 20 development samples | `results/` (`04_counts/`, `06_dge/`, `07_pathway/`, ...); 19 samples after QC (see QUANT step Q9) |
| Re-analyses, regression and validation runs | `results_rc1prep_19/` (git-ignored) |
| Reference (GENCODE v50, STAR/Salmon indexes) | `ref/` (git-ignored) |
| Public counts-level data (TCGA-COAD, GSE78220) | `data_public/` (git-ignored; URLs, sizes and sha256 in `data_public/SOURCES.tsv`) |
| Plan and pre-registered criteria | `results_rc1prep_19/pipeline_plan_20261015/PLAN.md` |
| Run history | `docs/run_log.md` (last entry 2026-10-15; the 2026-10-16 TCGA-COAD run is described in `results_rc1prep_19/validation_tcga_coad_diagnostic_20261016/REPORT.md` and is **not yet in the run log**) |

## 1. Entries

| Entry | Input | What it does | Real data |
|--|--|--|--|
| `-entry QUANT` | FASTQ sample sheet (csv) | reference, fastp, Salmon, STAR, featureCounts, per-sample QC, cohort QC, counts matrices, contract files, manifest | **NO (stub only)** in refactored form; the pre-refactor version ran on 20 samples |
| `-entry ANALYSIS` | counts matrix + sample sheet | input checks, DESeq2, ORA, GSEA, manifest | **YES** (19 samples, 2026-10-15) |
| `-entry DISCOVER --question diagnostic` | counts matrix + clinical file [+ hold-out list] [+ second route] | input contract, DESeq2, candidate rules 1-3b, LASSO score model, packaging | **YES** (19 samples, both routes; TCGA-COAD, one route) |
| SCORE (scripts, not Nextflow) | counts matrix | applies the locked model (`scripts/apply_locked_model.R`) and rule 4 (`scripts/rule4_reproduction.py`) | **YES** on the 19 development samples only; **NOT run on the 15 held-out samples** |
| default (no `-entry`) | FASTQ sample sheet | QUANT, then ANALYSIS fed by the contract files (as before the refactor; **not** DISCOVER) | NO (stub only) |

`--question prognostic` and `--question predictive` are **planned, not built**: DISCOVER stops with an error (`nextflow/main.nf` DISCOVER workflow; `config/questions.yaml:7-8`).

## 2. Step table

Parameter locations: `T:n` = `conf/thresholds.config` line n; `N:n` = `nextflow/nextflow.config` line n; `R:n` = `config/candidate_rules.yaml` line n; `Q:n` = `config/questions.yaml` line n. Other numbers are fixed in the named script.

### 2.1 QUANT (FASTQ to counts and QC)

| ID | Step | Input | Tool and version | Parameters and cutoffs (location) | Output | Check that proves it | Real data | Status |
|--|--|--|--|--|--|--|--|--|
| Q0 | Sample-sheet parse (`parseSamplesheet`, `main.nf:30-62`) | csv `sample,run,study,group` [`,batch,baseline`] | Nextflow 24.10.5 | `legacy_samplesheet = true` (N:26): 4-column sheet accepted, batch = study, baseline = group `normal`; run must match `^[SED]RR[0-9]+$`; exactly one baseline group | resolved sample rows | Stops with a named error on duplicates, empty group, bad accession, 0 or 2 baseline groups (code, `main.nf:47-59`). The legacy 4-column sheet ran for real; the refactored parser ran in stub tests A-F (`docs/run_log.md`, 2026-10-13). | YES (legacy sheet) | VERIFIED |
| Q1 | Memory guard (`resolveMethod`, `main.nf:75-83`) | cgroup memory limit | Nextflow | `min_mem_star_gb = 32` (T:16); `run_star = true` (N:17); `allow_salmon_only = false` (N:25) | method = STAR+featureCounts or Salmon-only, written to manifests | Below 32 GB (or unreadable) it stops with an error unless `--allow_salmon_only true`. Shown on this project's pods: stub test B stopped at 16 GB; `run_validation.sh preflight` FAILED on 32 GB (`docs/run_log.md`, 2026-10-13 and 2026-10-15). | YES (error path) | VERIFIED |
| Q2 | Reference build: `REF_FETCH`, `REF_FILTER`, `REF_GENOME_FILTER`, `REF_TABLES`, `SALMON_INDEX`, `STAR_INDEX` (`modules/reference.nf`, `modules/star.nf`) | GRCh38 FASTA, GENCODE GTF + transcripts | GENCODE v50 (`gencode_release`, N:13); salmon 1.10.3 (decoy-aware index, k = 31); STAR 2.7.11b | PAR1/PAR2 coordinates N:9-12; `star_index_ram = 48000000000` (N:22); `sjdbOverhang` = read length - 1 from fastp JSON | `ref/` incl. `gene_id_to_symbol.tsv` (78,670 genes) | GENCODE md5 checked against release MD5SUMS (`docs/tool_versions.tsv`). Built for real 2026-10-08/09. A `diff` against `modules/reference.nf.bak_pre_refactor` shows the refactor changed only the `stub:` text of `REF_TABLES` (line 122), so the real code is unchanged (VERIFIED by diff). | YES (pre-refactor; real code unchanged) | VERIFIED |
| Q3 | `VERIFY_FASTQ` | R1/R2 `.fastq.gz` | gzip test | `wait_hours = 48` (N:16); FASTQ path convention `data/<normal_sample or tumor_sample>/<run>_R1.fastq.gz` (`main.nf:147-150`) | readiness-gated reads | 20 of 20 samples passed for real (run log 2026-10-09 09:28). | YES (pre-refactor) | VERIFIED |
| Q4 | `FASTP` | trimmed reads | fastp 1.4.0 | `-q 20` (T:6), `--unqualified_percent_limit 40` (T:7), `--cut_tail` window 4 / mean q 20 (T:8-9), `--length_required 36` (T:10), `--n_base_limit 5` (T:11), `--detect_adapter_for_pe` (`modules/fastp.nf:23`, no number) | trimmed FASTQ, fastp JSON/HTML | 20 of 20 done, 0 failed (run log 2026-10-09). | YES (pre-refactor) | VERIFIED |
| Q5 | `SALMON_QUANT` | trimmed reads | salmon 1.10.3 | `-l A --validateMappings --gcBias --seqBias` (`modules/salmon.nf:18-19`; no numeric cutoff); executor memory 28 GB, `maxForks 2` (`nextflow.config:76`) | `quant.sf`, library type | 20 of 20 done (run log). Salmon peaks 25.9 / 26.6 GB (run log 2026-10-09). | YES (pre-refactor) | VERIFIED |
| Q6 | `STAR_ALIGN` | trimmed reads | STAR 2.7.11b | `--outFilterMultimapNmax 20` (T:14); one pass, unsorted BAM, other options STAR defaults (`modules/star.nf:53-57`); memory 33 GB, `maxForks 1` (`nextflow.config:83`) | BAM, `Log.final.out` | 19 of 19 started tasks done at 09:03; FEATURECOUNTS 20 of 20 at 09:28, which needs all 20 BAMs (run log). | YES (pre-refactor) | VERIFIED |
| Q7 | `SAMTOOLS_SORT_INDEX`, `SAMTOOLS_STATS` | BAM | samtools 1.24 | none numeric | sorted BAM, stats | 19 of 19 started tasks done at 09:03; FEATURECOUNTS 20 of 20 at 09:28 needs the sorted BAMs (run log). | YES (pre-refactor) | VERIFIED |
| Q8 | `FEATURECOUNTS` | sorted BAM + GTF | subread / featureCounts 2.1.1 | `-Q 10` (T:15); `-p --countReadPairs -B -C -t exon -g gene_id`; `-s` taken from the Salmon library type (`modules/featurecounts.nf:19-28`) | per-sample counts + summary | 20 of 20 done (run log). Strandedness consistent within study (study A reverse, study B unstranded; run log, step1_qc). | YES (pre-refactor) | VERIFIED |
| Q9 | `SAMPLE_QC_FLAGS` (`nextflow/bin/sample_qc.py`) | STAR log, featureCounts summary, Salmon dir, fastp JSON | python 3.12.15, pandas 2.3.3 | WARN: STAR unique < 70 % (T:19), featureCounts assigned < 60 % (T:20), Salmon mapped < 70 % (T:21); FAIL (sample excluded): STAR unique < 50 % (T:22), Salmon mapped < 50 % (T:23); new `reason` column (refactor) | `sample_qc.tsv` | Sample T5_B (SRR13401190) FAILED and was excluded, leaving the 19 development samples (run log 2026-10-09, step1_qc). The `reason` column was added after that run and was tested only with stubs. | YES (flags, pre-refactor); `reason` column NO (stub only) | VERIFIED (flags); NOT VERIFIED (`reason` column on real data) |
| Q10 | `MERGE_FEATURECOUNTS`, `SALMON_TXIMPORT` | per-sample outputs | python 3.12.15; tximport 1.38.2 (R 4.5.3) | tx2gene from the GENCODE v50 GTF | `counts_featurecounts.tsv`, `counts_salmon_tximport.tsv`, `salmon_tximport.rds` | Files exist in `results/04_counts/` and are the inputs of every later analysis (VERIFIED by use). `diff` against `modules/merge_counts.nf.bak_pre_refactor` shows the refactor changed only the `stub:` blocks (VERIFIED by diff). | YES (pre-refactor; real code unchanged) | VERIFIED |
| Q11 | `COHORT_QC` (`nextflow/bin/cohort_qc.R`) | counts, sample QC | R 4.5.3, DESeq2 1.50.2 | pre-filter >= 10 counts (T:26) in >= min(3, n) samples (T:27); outlier flag: robust z (MAD) > 3 on PC1 or PC2 within study x group cell (T:28), **flag only** | PCA, correlation heatmap, flags | 1 task done (run log 2026-10-09). The refactor only passes the three numbers from `conf/thresholds.config` to the script (`diff` against `modules/cohort_qc.nf.bak_pre_refactor`); the values are the old hard-coded ones (run log 2026-10-13, Step 3). | YES (pre-refactor) | VERIFIED |
| Q12 | `MULTIQC`, `QC_COPY_FASTP` | QC files | multiqc 1.35 | none | per-stage reports | 8 MultiQC tasks done (run log). | YES (pre-refactor) | VERIFIED |
| Q13 | `MAKE_CONTRACT`, `WRITE_MANIFEST` (`modules/contract.nf`, `bin/make_contract.py`, `bin/write_manifest.py`) | QUANT outputs | python 3.12.15 | `manifest_fastq_md5 = true` (N:28) | `contract/` (counts matrix, sample sheet, QC table), `manifest_quant.json` | `-stub-run` tests A, D, E, F passed; `write_manifest.py` run once for real on the reference and one FASTQ size (`docs/run_log.md`, 2026-10-13). **Never run on the real quantification.** | NO (stub only) | NOT VERIFIED on real data |

**QUANT as a whole:** the real quantification in `results/` came from the pre-refactor code. The refactored `main.nf` is NOT VERIFIED on real FASTQ: the 15-sample QUANT is blocked by the memory guard on the current pod.

### 2.2 ANALYSIS (counts to differential expression and pathways)

| ID | Step | Input | Tool and version | Parameters and cutoffs (location) | Output | Check that proves it | Real data | Status |
|--|--|--|--|--|--|--|--|--|
| A1 | `VALIDATE_ANALYSIS_INPUT` (`bin/validate_analysis_input.py`) | counts matrix, sample sheet, symbols, optional QC table | python 3.12.15 | >= 3 samples per group (T:31, PROVISIONAL, no source); <= 5 % gene IDs unmapped (T:32, PROVISIONAL, no source); integer counts only when `count_source = featurecounts`; one baseline group; group not confounded with batch | `sample_sheet.validated.tsv`, JSON of checks | Negative tests give the intended `CHECK` line and exit 2 (`docs/run_log.md`, 2026-10-13); passed on the 19 samples (2026-10-15). The two thresholds are PROVISIONAL: the numbers have no external source. | YES | VERIFIED (code path); the thresholds themselves are NOT VERIFIED (no source) |
| A2 | `DESEQ2` (`bin/deseq2.R`) | counts, validated sheet | R 4.5.3, DESeq2 1.50.2, apeglm 1.32.0 | design `~ batch + group` (`~ group` if one batch), reference = baseline group from the sheet (`deseq2.R:14-19`); pre-filter >= 10 counts (T:35) in >= n samples, n = `'auto'` = smallest group (T:36); `results(alpha = 0.05)` (T:37, same value as the padj cut); significance table uses padj < 0.05 and **unshrunk** \|log2FC\| >= 1 (T:38; shrunk-rule switch not done); `lfcShrink(type = "apeglm")` (`deseq2.R:27`); 15 volcano labels (T:39), 50 heatmap genes (T:40) | `deseq2_results_all.tsv`, significant up/down tables, MA, volcano, heatmap, summary JSON | Re-run via `-entry ANALYSIS` on the 19 samples: 5 of 5 tasks, 9 min 50 s, results **byte-identical** to the earlier fit (`docs/run_log.md`, 2026-10-15 b). | YES | VERIFIED |
| A3 | `PATHWAY_ORA` (`bin/ora.R`) | DESeq2 table | clusterProfiler 4.18.4, org.Hs.eg.db 3.22.0, msigdbr 26.1.1 | significant set: padj < 0.05, unshrunk \|log2FC\| >= 1 (`de_padj`, `de_lfc`); universe = tested genes with a symbol; set size 10-500 (T:43-44); BH q 0.05 (T:45); GO BP; KEGG **off** by default (`run_kegg = false`, N:39; licence unresolved); 15 plot categories (T:48) | `ora_*` tables and dotplots | 12 DESeq2/ORA/GSEA tables identical between the one-command run and run b (`docs/run_log.md`, 2026-10-15 f). Also run on the TCGA-COAD fit (exit 0). | YES | VERIFIED |
| A4 | `PATHWAY_GSEA` (`bin/gsea.R`) | DESeq2 table | fgsea 1.36.2, BiocParallel 1.44.0, msigdbr 26.1.1 | ranking = Wald stat, one row per symbol; `fgseaMultilevel`, `nPermSimple 10000` (T:46), seed 42 (T:47), padj 0.05 (T:45), set size 10-500; collections Hallmark and GO BP (KEGG only with `run_kegg true`); 20 plot sets (T:49) | `gsea_*` tables and plots | as A3. **With KEGG off, GO BP draws a different position of the random stream, so its permutation-based values can differ slightly from a run with KEGG** (comment at the top of `gsea.R`). | YES | VERIFIED |
| A5 | `SECONDARY_ROUTE`: Salmon DESeq2 + `DE_CONCORDANCE` (default entry only) | Salmon tximport | as A2; `concordance.R` | same `de_padj`, `de_lfc` | `deseq2_salmon_*`, concordance tables | Ran in the real pre-refactor run (09:28, 1 task each). Not part of `-entry ANALYSIS`. | YES (pre-refactor) | VERIFIED |
| A6 | `WRITE_MANIFEST` (analysis) | all inputs and outputs | python 3.12.15 | 60 resolved parameters written | `manifest_analysis.json` | `results_rc1prep_19/diagnostic_final_20261015/analysis_entry_19/out/contract/manifest_analysis.json` exists (VERIFIED by `find`; contents not re-read) | YES | VERIFIED (file exists) |

### 2.3 DISCOVER (question-driven discovery; diagnostic only)

Entry: `nextflow/main.nf` workflow `DISCOVER`, processes in `modules/discover.nf`. It never writes into `lock/`.

| ID | Step | Input | Tool and version | Parameters and cutoffs (location) | Output | Check that proves it | Real data | Status |
|--|--|--|--|--|--|--|--|--|
| D1 | `DISCOVER_INPUT_CHECK` (`bin/check_discover_input.py`) | counts, clinical file, symbols, optional hold-out and second route | python 3.12.15, pandas 2.3.3 | C0 question implemented; C1 required clinical columns; C2 sample IDs one-to-one; C3 raw integer counts; C4 two groups and one baseline; C5 hold-out is patient-level (no patient in both sets); C6 >= 3 samples per group (T:31, PROVISIONAL); C7 group not confounded with batch; C8 <= 5 % unmapped (T:32, PROVISIONAL); C9 leave-one-out only up to 60 samples (`Q:29`, PROVISIONAL); warnings W1-W4 (no batch, repeated patients, single batch, single route) | `input_checks.json`, training sample sheet and counts | 10 negative tests each stop at the intended check, valid input passes (`discover_regression_20261015/input_check_tests/`, 2026-10-16 run-log entry). TCGA-COAD run: PASS with 4 warnings. | YES | VERIFIED |
| D2 | `DISCOVER_DESEQ2` (one or two routes; same `deseq2.R` as A2) | training counts | as A2 | same parameters as A2 | `01_de/` tables | 19-sample regression: DESeq2 tables of both routes byte-identical to the earlier fits; TCGA re-run byte-identical (D5). | YES | VERIFIED |
| D3 | `CANDIDATE_RULES` (`bin/candidate_rules.R`) | counts, DESeq2 tables, `config/candidate_rules.yaml`, `config/questions.yaml` | R 4.5.3, DESeq2 1.50.2, apeglm 1.32.0 | rule 1: padj < 0.05 (R:12) and \|shrunk log2FC\| >= 1 (R:13); rule 2: baseMean >= 20 (R:17; stricter flag 50, R:19-20, report only); rule 3a: rule 1 on both routes, same direction (R:22-25); rule 3b: `loo` (`Q:27`, locked method for 19 samples; pass >= `ceil(0.90 x refits)` = 18 of 19, R:32, `Q:28`) or `subsample` (20 patient-level 80 % refits, `Q:30-31`, seed `Q:32`; PROVISIONAL). **Single route:** rule 2 applied to the one route, rule 3a SKIPPED (`Q:19-20`, `candidate_rules.R:102-108`), and the report says candidates are not equivalent to the 19-sample method | `candidate_table.tsv`, funnel, pass matrix | 19-sample regression: candidate table, pass matrix and funnel byte-identical to the locked run; 1,390 candidates (510 up, 880 down, 1,353 protein-coding). The code stops if the derived pass count differs from the locked 18. TCGA-COAD (single route, subsample): 3,775 candidates. | YES | VERIFIED |
| D4 | `FIT_MODEL` (`bin/fit_diagnostic_model.R`) | counts, candidates | R 4.5.3, glmnet 4.1.10, pROC 1.19.1 | features = candidates passing all applied rules; log2(CPM + 1), CPM total = column sum (`fit_diagnostic_model.R:21`); `cv.glmnet` binomial, alpha 1, deviance, `lambda.1se` (line 32), K = 5 (`Q:37`) folds stratified by group x batch, seed 42 (`Q:35`); cutoff = median of 5 inner-fold Youden cutoffs, seed 43 (`Q:36`) (lines 36-37); score = mean up genes - mean down genes; call tumour if score > cutoff | `model_genes.tsv`, `model_spec.txt`, development scores (training, optimistic) | 19-sample regression: 9 genes, signs, coefficients, cutoff 1.352325902 and scores identical to the lock (difference 0). | YES | VERIFIED |
| D5 | `PACKAGE_MODEL` (`bin/package_model.py`) | model files, checks | python 3.12.15 | none | `MODEL_MANIFEST.tsv`, `DISCOVER_SUMMARY.md` (labels results "development, optimistic; not the locked model") | produced in all DISCOVER runs | YES | VERIFIED |
| D6 | Hold-out scoring, pathways, comparison with a published list | model files, hold-out list, DE table | R scripts in the validation folder (`score_holdout.R`, `compare_hosseini.py`), `ora.R`/`gsea.R` | hold-out scored once; CIs: DeLong (AUC), Clopper-Pearson (sens., spec.) | metrics tables | Done for TCGA-COAD (`validation_tcga_coad_diagnostic_20261016/`). **MANUAL: not a step of the Nextflow DISCOVER entry**, so hold-out scoring and pathways must be run by hand. | MANUAL | VERIFIED (TCGA run); NOT in the entry |

### 2.4 SCORE (apply the locked model; scripts)

| ID | Step | Input | Tool and version | Parameters and cutoffs (location) | Output | Check that proves it | Real data | Status |
|--|--|--|--|--|--|--|--|--|
| S1 | `scripts/apply_locked_model.R` | featureCounts matrix (versioned GENCODE v50 IDs), optional sheet and QC table | R 4.5.3, pROC 1.19.1 | reads genes, signs and cutoff 1.352325902 from `lock/diagnostic_v1/` after checking `sha256sums_final.txt`; stops if a model gene is missing; QC FAIL samples dropped and listed; `--evaluate` on validation data requires the QC table | scores, calls, optional performance | Scores equal `lock_proposal_v2/development_scores.tsv` (difference 0) on the 19 development samples (2026-10-15 c); 5 negative tests stop with exit 2. **Development data only, optimistic.** Not run on any independent sample. | YES (development only) | VERIFIED (development); NOT VERIFIED (independent data) |
| S2 | `scripts/rule4_reproduction.py` | DESeq2 table | python 3.12.15 | same direction and padj < 0.05 (`de_padj`, T:37); no pass threshold (R:34-39) | share of the 1,390 candidates reproduced | Self-consistency on the 19 samples: 1,390 of 1,390 (both routes). **Not run on independent samples.** | YES (self-consistency only) | VERIFIED (self-consistency); NOT VERIFIED (independent) |
| S3 | `scripts/run_validation.sh` (`preflight`, `quant`, `resume`, `downstream`) | 15 held-out samples (`validation_run/samplesheets/validation_15.csv`) | Nextflow, S1, S2 | preflight needs cgroup memory >= 60 GB and >= 16 cores; pass rule at most 1 error per class, overall and per study (R:45) | validation QC, locked-model calls, rule 4 | Preflight FAILED on the 32 GB / 8-core pod (2026-10-15) and would fail on this 16 GB / 4-core pod. 30 of 30 FASTQ present with ENA sizes (96.1 GB). **QUANT, scoring and rule 4 were never run.** | NO | NOT VERIFIED |

### 2.5 Manual steps and regression tests (outside the entries)

| Script | Purpose | Status |
|--|--|--|
| `scripts/candidate_selection.R` | the original manual rules 1-3b (DISCOVER's `candidate_rules.R` is its generalised copy); guarded to the 19 samples | MANUAL; ran; DISCOVER reproduces it exactly |
| `scripts/downstream/reactome_analysis.R`, `reference_gene_check.R`, `lasso_stage1_deseq2.R`, `lasso_stage2_glmnet.R`, `lasso_arm3.R`, `lasso_arm4.R`, `relative_expression_dev.R` | Reactome ORA/GSEA, reference genes, LASSO trial arms. Some thresholds are hard-coded copies of the config (e.g. reference-gene baseMean 100, `reference_gene_check.R:13`) | MANUAL; tied to the 19 samples; **not** part of any entry |
| `results_rc1prep_19/lock_proposal_v2/fit_final_model_v2.R` | the original LASSO fit that produced the lock | MANUAL; DISCOVER's `fit_diagnostic_model.R` reproduces it |
| `scripts/run_diagnostic.sh` | one command: ANALYSIS, locked model, rule 4, regression | ran 2026-10-15 (541 s, 48 of 48 PASS) |
| `tests/regression/check_19_development.py`, `check_discover_regression.py` | 48 expected-value checks on the 19 samples; DISCOVER-vs-lock checks | 48 of 48 PASS; DISCOVER 23 of 23 IDENTICAL or PASS |
| `validation/check_vs_paper.*` | comparison with Hosseini 2023 (validation only, not part of QUANT or ANALYSIS) | MANUAL |

## 3. How to run

Run from `/workspace/Biomarker_pipeline`. Use a new `--outdir` and a new work directory for every run. The wrapper scripts (`run_diagnostic.sh`, `run_validation.sh`, `apply_locked_model.R`) refuse an existing output folder; plain `nextflow run` does not refuse one, so overwrite protection there is by convention (NOT VERIFIED that DISCOVER refuses an existing `--outdir`). Nextflow is `env/bin/nextflow` (24.10.5).

**3.1 FASTQ input (QUANT).** Sample sheet columns: see 4.1. FASTQ files go in `data/normal_sample/` (baseline group) and `data/tumor_sample/` (other group) as `<run>_R1.fastq.gz`, `<run>_R2.fastq.gz`.
```
env/bin/nextflow run nextflow/main.nf -entry QUANT -profile conda --samplesheet <sheet.csv> --outdir <new dir> -w <new work dir> -with-trace <trace.tsv>
```
Output: counts matrices, QC tables and `contract/` (`counts_matrix.tsv`, `sample_sheet.tsv`, `qc_table.tsv`, `manifest_quant.json`). For the 15 held-out samples use `scripts/run_validation.sh preflight` then `quant` (the memory guard is never bypassed). **NOT VERIFIED on real FASTQ in the refactored form.**

**3.2 Counts input, differential expression and pathways (ANALYSIS).**
```
env/bin/nextflow run nextflow/main.nf -entry ANALYSIS -profile conda --counts_matrix <counts.tsv> --sample_sheet <sheet.tsv> --count_source featurecounts --outdir <new dir> -w <new work dir>
```
Add `--count_source salmon_tximport --counts_rds <salmon_tximport.rds>` for the Salmon route, `--qc_table`, `--quant_manifest` as needed. One command for the 19-sample development set: `scripts/run_diagnostic.sh results_rc1prep_19/diagnostic_run_<YYYYMMDD_HHMM>`.

**3.3 Counts input, discovery (DISCOVER).**
```
env/bin/nextflow run nextflow/main.nf -entry DISCOVER -profile conda --question diagnostic --counts_matrix <featureCounts counts.tsv> --clinical <clinical.tsv> \
    [--holdout <sample ids, one per line>] [--counts_matrix_secondary <salmon counts.tsv> --counts_rds_secondary <tximport.rds>] \
    [--questions_config <copy of config/questions.yaml>] --outdir <new dir> -w <new work dir>
```
- With one route only, rules 2 and 3a are modified or skipped (D3). For more than 60 training samples, copy `config/questions.yaml`, set `robustness_method: subsample`, and pass the copy with `--questions_config` (the repository copy is not edited).
- A hold-out with two routes is **not supported** (check C5 fails).
- On a pod smaller than the default executor (16 cores, 58 GB: `nextflow.config:61`), add a `-c` config that lowers the executor and process memory (as done in `validation_tcga_coad_diagnostic_20261016/run_config/pod4.config`).

**3.4 FASTQ input to DISCOVER.** There is no single command: run QUANT (3.1), then DISCOVER (3.3) on `counts_featurecounts.tsv` plus a clinical file you write (columns in 4.2). The default entry (no `-entry`) runs QUANT then ANALYSIS, not DISCOVER. NOT VERIFIED end to end.

**3.5 Score with the locked model.**
```
env/prefix/r_classifier/bin/Rscript scripts/apply_locked_model.R --counts <counts.tsv> --cohort development|validation [--sample-sheet <sheet.tsv>] [--qc-table <qc.tsv>] [--evaluate] --out-dir <new dir>
env/prefix/pyqc/bin/python3 -I scripts/rule4_reproduction.py ...   # see the script header for arguments
```
The locked model is for featureCounts matrices on **versioned GENCODE v50 gene IDs** (data from other annotations need ID mapping first).

## 4. Required input files

**4.1 QUANT sample sheet (csv)**

| Column | Required | Rule |
|--|--|--|
| `sample` | yes | unique label |
| `run` | yes | SRA accession `^[SED]RR[0-9]+$` |
| `study` | yes | study accession |
| `group` | yes | any two group names (not fixed in code); exactly one baseline group |
| `batch` | no | design term; default = study (only with `legacy_samplesheet = true`, N:26) |
| `baseline` | no | `TRUE` / `FALSE`, one baseline group; default = group named `normal` (legacy only) |

**4.2 DISCOVER clinical file (TSV, one row per sample)**

| Column | Required | Rule |
|--|--|--|
| `sample_id` | yes | equals the matrix column names, unique |
| `patient_id` | yes | used for patient-level hold-out and flagged when repeated |
| `group` | yes | exactly two levels |
| `baseline` | yes | `TRUE` for exactly one group |
| `batch`, `study`, `label` | no | `batch` (or `study`) becomes a design term; none = one batch, design `~ group` |

**4.3 ANALYSIS sample sheet (contract format, TSV)**: `sample_id, label, group, baseline, study, batch, count_source`. **Counts matrix:** genes (rows, Ensembl IDs) x samples (columns); raw integer counts for `featurecounts`, tximport counts for `salmon_tximport`. **Reference table** `ref/gene_id_to_symbol.tsv`: `gene_id, gene_name, gene_type` (GENCODE v50).

## 5. Resources

| Use | Requirement | Source |
|--|--|--|
| QUANT (STAR route) | cgroup memory >= 32 GB for the guard (T:16); in practice >= 60 GB and 16 cores (STAR ~28.5-33 GB x1 plus Salmon 28 GB x2 inside a 58 GB executor budget, `nextflow.config:61,76,83`); FASTQ ~44 GB for 20 samples, 96 GB for the 15 held-out | `run_validation.sh:18-19`; run log 2026-10-09 (VERIFIED on the 64 GB / 16-core pod) |
| ANALYSIS (19 samples) | < 2 GB, about 10 min on 8 cores | `docs/README_diagnostic_workflow.md`; run 541 s (VERIFIED) |
| DISCOVER (19 samples, both routes, leave-one-out) | 24 min on this project's 8-core pod | run log entry 2026-10-16 (VERIFIED) |
| DISCOVER (TCGA-COAD, 332 training samples, subsample) | about 23-24 min per run on 4 cores / 16 GB (two runs) | `validation_tcga_coad_diagnostic_20261016/run1/timing.txt`, `run2/timing.txt` (VERIFIED) |
| Software | micromamba envs in `env/prefix/` built from `env/*.yml`; versions in `docs/tool_versions.tsv`; GENCODE v50 | VERIFIED |

## 6. Limitations and open items

1. **Not clinical validation.** Everything here is technical or development work. No clinical claim is made or supported.
2. **QUANT on real FASTQ was not re-run after the refactor.** The real quantification (20 samples) came from the pre-refactor code; the refactored QUANT, contract and manifest steps have only been stub-tested. NOT VERIFIED on real data.
3. **The 15-sample independent validation has not been run.** The pod is too small (needs >= 60 GB and 16 cores). The locked 9-gene model and rule 4 therefore have no independent validation; all locked-model numbers are development, optimistic.
4. **Unstranded assumption (TCGA-COAD).** The `unstranded` GDC count column was chosen on the assumption that the libraries are unstranded. NOT VERIFIED.
5. **Single-route rules skipped.** With one quantification route, rule 3a is skipped and rule 2 is applied to one route, so candidates are not equivalent to the 19-sample method. Rule 3b on large cohorts uses PROVISIONAL subsampling values (`Q:30-32`) instead of leave-one-out.
6. **No patient or batch term in the TCGA-COAD run.** The design was `~ group`; 27 patients have a tumour and a normal sample in training and the model has no patient term. Site or plate effects were not modelled.
7. **Commercial-use terms for TCGA (GDC) and GEO data are open.** GDC pages do not address commercial use, and the GEO terms page could not be read (`PLAN.md` section c). NOT VERIFIED. KEGG is off by default for a similar licence reason.
8. **Prognostic and predictive modes are planned, not built.** `--question prognostic|predictive` stops with an error; `config/questions.yaml:40-41` holds PROVISIONAL placeholders only.
9. **PROVISIONAL thresholds without a source:** `analysis_min_samples_per_group = 3`, `analysis_max_unmapped_gene_pct = 5` (T:31-32); `robustness_loo_max_samples = 60`, `robustness_subsamples = 20`, `robustness_subsample_fraction = 0.80` (`Q:29-31`).
10. **Known gaps in the tooling:** hold-out scoring and pathways are not Nextflow steps of DISCOVER; no container, CI or git remote; the DESeq2 significance table still uses the unshrunk log2FC (the candidate rules use the shrunk value); `docs/run_log.md` lacks the 2026-10-16 TCGA-COAD entry; the owner's approval text for `PLAN.md` was not recorded on disk.
