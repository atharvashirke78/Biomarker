# Upgrade backlog: diagnostic RNA-seq biomarker pipeline

Written 2026-10-22 from the validation audit (`/workspace/run_dev19/validation_audit_20261022T014443Z/`, findings `AUD-*` in `findings.tsv`, code review `CR-F01…F23`) and the work done since. **Nothing in the "Open" sections has been applied.** Items marked **[scientific]** change an analysis method, threshold or selection rule and need the owner's approval first. Items marked **[technical]** change behaviour around the analysis (checks, messages, structure) and must still be re-tested before use. The real-data validation run of the standardized scripts (19 samples, output `result_final`) was in progress when this was written; reprioritise after its comparison with the earlier run.

Status key: **Done** (applied this session), **Open** (not started), **Decision** (needs an owner decision before any work).

## A. Done this session (for the record)
| ID | What | Files |
|---|---|---|
| D1 | Shipped config has no dev paths in the fixed block; input/output paths are the two user settings at the top (now set for the validation run) and still overridable with `--input_dir` / `--output_dir` (tested: derived paths follow the override) | `nextflow/nextflow.config` |
| D2 | Mate-tag styles `_R1_001/_R2_001`, `_R1/_R2`, `.R1/.R2`, `_1/_2`; group aliases case/control, disease/healthy added | `conf/input.config` |
| D3 | featureCounts stops when the Salmon library format cannot be read (no silent `-s 0`); QC merge uses `pipefail` and `LC_ALL=C` | `modules/featurecounts.nf`, `modules/sample_qc.nf` |
| D4 | The 19 leave-one-out refits get the same DESeq2 settings as the main fit (resolved values passed in) | `nextflow/bin/candidate_rules.R`, `modules/discover.nf` |
| D5 | Manifest records real tool versions (was "NOT AVAILABLE") | `nextflow/bin/write_manifest.py`, `modules/contract.nf` |
| D6 | CPU and RAM budget follows the pod; parallel STAR and Salmon jobs derived from it (same values as before on the 64 GB pod) | `nextflow/nextflow.config` |
| D7 | Report and index text no longer say "colon"; the early input scan prints non-fatal, per-stage warnings for small (< 10 leave-one-out refits) and large (> 60, leave-one-out) cohorts; "NO MODEL FITTED" message explains rule 3b not evaluable | `scripts/make_report.py`, `nextflow/main.nf`, `nextflow/bin/package_model.py`, `scripts/input_scan.py`, `modules/folder_input.nf`, `nextflow/bin/fit_diagnostic_model.R` |
| D8 | Documentation corrected against the code (required vs optional metadata columns, labels, naming, `--check_only`, output paths); `docs/INTEGRATION_GUIDE.md`; README status and run times | `README.md`, `docs/*.md` |
| D9 | Thin job service (upload, validate with the pipeline's own scanner plus a service check for unsafe batch/study values, one-at-a-time queue, live workflow chart, summary, figures, resume, delete, report download); no token or login (localhost only, remote bind refused unless `SERVICE_ALLOW_REMOTE=1`); demo button with real data in `/workspace/biomarker_demo`, live animated workflow view with step descriptions, About page; 21 tests | `service/` |

## B. Before the pipeline can be called a release candidate
| ID | Upgrade | Type | Why (evidence) | Files / where | Verify by |
|---|---|---|---|---|---|
| U1 | Commit the reviewed code on `feat/one-workflow`, tag it, regenerate the freeze list, complete the backup of changed files; record the commit in the freeze record | technical, needs your go-ahead | Code ran uncommitted; run_info shows commit `46c4860` that does not identify the code (AUD-02); `docs/FROZEN_FILES.standardized.sha256` is stale | repo root, `docs/FREEZE_RECORD_DRAFT.md` | clean `git status`, hash list verifies |
| U2 | Compare the validation run with the earlier run against pre-agreed criteria (count matrix hash `cce23f518d50…`, featureCounts-route DESeq2 table, 19 refits, cutoff 1.352325902; allowed differences on the Salmon route for the four known outlier samples) | technical | Exit code 0 is not validation | `/workspace/run_dev19/compare_with_october.py`, `step_compare.py` | signed comparison workbook |
| U3 | Report honesty fixes: write "completed" only after the run ends (AUD-05); make the "MECHANICS TEST ONLY" limit a documented parameter; state unshrunk vs shrunk log2FC where each is used; add sections "limitations", "uncertainty", "recommended validation steps"; keep generated text separate from interpretation | technical | `make_report.py`, report banner shown for a 9-tumour group; two LFC definitions (below) | `scripts/make_report.py`, `main.nf` onComplete | read the final report |
| U4 | One definition of the DE log2FC cutoff: pathway/volcano/significant tables use the unshrunk value (`deseq2.R:34-38`, `ora.R:11`, `de_lfc`), candidate rule 1 uses the apeglm-shrunk value (`candidate_rules.yaml`) | **scientific** (choice of quantity) | Different gene lists feed pathways and candidates | `conf/thresholds.config`, `nextflow/bin/deseq2.R`, `ora.R`, `config/candidate_rules.yaml` | counts of significant genes before/after |
| U5 | Document and mitigate Salmon run-to-run bimodality: record it in the report; options: rule 3a margin, repeated Salmon runs with a consensus, or decide on the featureCounts route only | **scientific** | AUD-11: 4 of 5 re-run samples fall into two solutions (about 2 % of reads); 19 rule-1 flips, 14 candidate-table changes, one final candidate (CDK2) at risk. **Single-thread test (2026-10-22, same trimmed reads of SRR10043499): two `-p 1` runs still differ by 2.6 % of reads at transcript level (two `-p 8` runs: 4.0 %), and Salmon 1.10.3 shows no seed option, so thread count is not the cause and `-p 1` does not fix it** | `config/candidate_rules.yaml`, `candidate_rules.R`, Salmon module | repeat-run experiment on all 19 samples |
| U6 | Model stability: fixed (sorted) feature order, stability selection or repeated cross-validation, report selection frequencies; nested cross-validation for an honest performance estimate | **scientific** | AUD-12: alphabetical order selected 8 genes (ETV4 lost); the lock records 21 of 30 seeds giving the same 9 genes; development performance is in-sample | `nextflow/bin/fit_diagnostic_model.R`, `config/questions.yaml` | selection frequency table, nested-CV AUC with interval |
| U7 | Independent validation: score the 15 held-out samples once with the locked model; hold-out support in the one-workflow route; plan for public cohort GSE142279 (40 paired runs; raw reads only, strandedness, site and purity unverified) | **scientific** | AUD-22, AUD-23: no independent validation exists | `scripts/apply_locked_model.R`, `docs/VALIDATION_LEDGER_DRAFT.md`, `docs/independent_validation_prompt.md` | pre-registered pass rule, one application |
| U8 | Add the scoring code and selection scripts to the lock hash set and correct the "DRAFT (not locked)" header in the spec | **scientific/owner** (changes the lock) | AUD-18: 2 of 7 lock files hashed at scoring time, no code hashed | `lock/diagnostic_v1/` | new hash list, apply test |

## C. Robustness of the workflow (technical, no change of results expected)
| ID | Upgrade | Evidence | Files |
|---|---|---|---|
| U9 | Gate the cohort steps on `CHECK_QC` so nothing runs on too few samples | CR-F05 | `nextflow/main.nf` |
| U10 | `NO MODEL FITTED` should be visible in the workflow status (banner, exit status option), not only in files | AUD-09 | `fit_diagnostic_model.R`, `main.nf`, `make_report.py` |
| U11 | `failOnMismatch` on per-sample joins; check that matrix columns equal the non-FAIL samples | CR-F02, F10, F11 | `main.nf`, `make_contract.py` |
| U12 | `--check_only` (or a new `--preflight`) should also run the input scan and not create the output folder or logs | AUD: check_only checks software and presence only, then `onComplete` creates `<output>/logs` | `main.nf`, `modules/folder_input.nf`, `nextflow.config` trace blocks |
| U13 | Graceful handling above 60 samples (clear skip of the discovery stage with a report) instead of failing after quantification | `check_discover_input.py` C9 | `check_discover_input.py`, `main.nf` — **Decision**: the alternative method (`subsample`) is a scientific choice |
| U14 | Duplicate gene-ID check in the locked scorer should use its own message and exit code | AUD-17 | `scripts/apply_locked_model.R` |
| U15 | `candidate_rules.R` should receive `de_padj` etc. only as resolved values (remove the text re-read) | CR-F06 (partly done in D4) | `candidate_rules.R` |
| U16 | Salmon memory scaling on retry, retry exit-code list | CR-F19 | `nextflow.config` |
| U17 | Launch from a dedicated launch folder (not the code checkout), keep `.nextflow*`, `run1.log` and `result/` out of git | AUD-24; the validation run was launched from `nextflow/` | `.gitignore`, run docs |
| U18 | Process-level isolation for long runs: start with `setsid nohup` or a service, never as a job of an interactive terminal (a closed terminal kills the run) | the resumed validation run had SIGHUP not ignored; an earlier attempt was killed by an outside signal | run docs, `service/` |
| U25 | **Bug found in the per-module check (2026-10-22):** a `batch_or_study` value with a space (for example `Study A`) passes the input scan, then `SAMPLE_QC_FLAGS` crashes with `ValueError` after the whole per-sample phase (hours), because `meta.study` is put unquoted into the shell command (`modules/sample_qc.nf:15`). A comma (`A, B`) is accepted by the scanner but stops the run early with "Bad samplesheet row" (`nextflow/main.nf:48`, naive comma split; `main.nf:158` writes the resolved sheet unquoted). Smallest fix: make the scanner reject batch/study values outside `[A-Za-z0-9._-]` with a clear message (same rule as sample IDs), and quote `${meta.study}` / `${meta.group}` in `sample_qc.nf` as a second guard | technical, input-validation rule: needs your approval | reproduced read-only with yesterday's real task inputs | `scripts/input_scan.py`, `modules/sample_qc.nf`, tests |
| U19 | Update stale docs (`docs/two_pipeline_design.md`, the `scripts/rnaseq.sh:8` comment, the legacy wrapper route) or retire `scripts/rnaseq.sh` | docs audit | `docs/`, `scripts/rnaseq.sh` |

## D. Broader input support (clinical and preclinical use)
| ID | Upgrade | Type | Notes |
|---|---|---|---|
| U20 | Merge multiple lanes per sample automatically; accept uncompressed `.fastq`/`.fq`; single-end input | technical | Today lanes become separate samples and other extensions are ignored |
| U21 | Other species (mouse, rat): parameterise reference, GTF, gene-ID mapping, pathway collections | technical + scientific | Pipeline is human GRCh38 / GENCODE v50 only; preclinical research needs this |
| U22 | More than two groups, covariates (age, sex, batch beyond one factor), paired design (patient term) | **scientific** | Design is `~ batch + group`; paired tumour/normal patients are not modelled as pairs |
| U23 | Support count-matrix input in the one-workflow route (the scanner already has `--counts`) | technical | For public datasets such as TCGA |
| U24 | Extra quantification options (salmon only route) with full validation | **scientific** | Salmon-only route exists in code but the one-workflow entry rejects it |

## E. Throughput and service
| ID | Upgrade | Notes |
|---|---|---|
| S1 | Concurrency by pod size is already derived (Salmon jobs = min(RAM / 28 GB, CPUs / 8)); test on a larger pod (for example 128 GB, 32 CPUs) | **NOT TESTED** on a bigger pod |
| S2 | Several pods, each with its own service and queue; shared storage and one job database | needed only for many cohorts |
| S3 | Cluster or cloud executor profile (Slurm, Kubernetes, AWS Batch) with the reference on every node | needs re-testing for identical results |
| S4 | Service hardening: user accounts and per-user isolation, TLS, retention and deletion rules, audit log, upload resume, job resume, progress percentage | thin version is local-only; optional shared token, Host/Origin guard and sandboxed reports added 2026-10-23; user accounts, TLS, audit log still open |
| S5 | Reduce STAR start-up cost per sample (shared index in memory) | measure first; STAR is 1 job at a time today |

## F. Methodological contribution (do not label as novel yet)
- Prognostic and predictive objectives need survival, event and treatment data that the current metadata and clinical file do not carry (`config/questions.yaml` marks both as not implemented; `docs/PREREG_prognostic_DRAFT.md` is a draft). Each needs its own pre-specified method, sample-size rule and validation design.
- The standardized cutoff procedure (median of inner-fold Youden cutoffs, fixed seeds) and the reproducible end-to-end design must be compared with published approaches before they are called a methodological contribution. No literature review has been done; reports and capability statements must not use "novel" until it has.
- Training, cross-validation and independent validation must stay separate: a candidate list from a development cohort is a hypothesis, not a validated biomarker.

## Suggested order
1. Finish and compare the validation run (U2), fix report wording (U3), commit and tag (U1).
2. Decide U4–U6 (scientific); run the independent validation once (U7).
3. Robustness items C (U9–U19), then service hardening (S4).
4. Broader inputs (D) and prognostic/predictive work only after the diagnostic pipeline is released.

## G. Check of demo job `5d6ede1e1e43` and comparison with the earlier run of the same samples (2026-10-22)
The demo job (6 real samples, 100,000 read pairs each; service demo button; 74 tasks, 10 min 34 s, exit 0) was compared with the earlier run of the identical inputs, `/workspace/example_discover_20261021` (earlier code), with the new read-only script `tests/regression/compare_runs.py` (result table: `validation_audit_20261022T014443Z/demo_vs_example_comparison.tsv`). Usage: `<pyqc python> tests/regression/compare_runs.py --new <folder> --base <folder>`; exit code 1 when a must-be-exact check differs.

| Check | Result |
|---|---|
| featureCounts count matrix (sha256) | **identical** (`b84c0add494b68d9`) |
| sample QC table | **all cells equal** |
| DESeq2, featureCounts route (06_dge and 09_discover/01_de) | **1236 genes x 9 columns, every cell equal** |
| leave-one-out pass matrix, featureCounts columns of the candidate table, candidate funnel after rule 1 on featureCounts | **identical** |
| manifest tool versions (fix D5) | **now recorded in a real run** (fastp 1.4.0, STAR 2.7.11b, Salmon 1.10.3, featureCounts v2.1.1, DESeq2 1.50.2, ...) |
| report text old vs new | only the intended changes (new warning, intended-use line, banner) |
| Salmon route | differs: gene counts 0.02-0.09 % of counts per sample (3-18 genes differ by > 1 read), but DESeq2 statistics change for every gene (max shrunk log2FC difference 1.25, baseMean up to 51, padj for 734 of 1424 genes); one gene (EHF, `ENSG00000135373.15`, padj 0.024, baseMean 14.3) dropped out of rule 1 on Salmon: funnel 20 to 19, candidate table 61 to 60 rows; no later stage changed |

Conclusion: the featureCounts-based path is exactly reproducible on this input; the only difference is Salmon's non-repeatability (U5). Biology of the figures is plausible (volcano: GUCA2A, CA4, AQP8, MYL9 down in tumour, DUOX2 up), but **PC1 (71 % of variance) follows the study, not the group**.

### Report and input defects seen in this job (all in `scripts/make_report.py` unless stated; not fixed)
| ID | Defect | Where | Fix (smallest) |
|---|---|---|---|
| R1 | Section "2. Stages" is empty in the one-workflow route | `make_report.py:153-158` (reads `nf_<stage>/exit_code.txt`, absent here) | list stages from the trace, or drop the heading |
| R2 | A stray `_(none)_` is printed after "NO MODEL WAS FITTED" | `:191` (model-genes table of an unfitted model) | print the table only when a model was fitted |
| R3 | Stale sentence "Pathway analysis inside DISCOVER (Reactome/GSEA port) is NOT built yet" although ORA and GSEA ran | `:194` (CR-F16) | delete or replace by the real pathway summary |
| R4 | Banner "Gene selection was done on all development samples ..." is shown when no selection happened | `:119`, `:127` | show it only when a model was fitted |
| R5 | "NO MODEL WAS FITTED: fewer than 2 candidate genes passed" does not say that rule 3b was not evaluable (the real reason) | `:190` | reuse the reason written by `fit_diagnostic_model.R` (D7) |
| R6 | Two DE sections with identical plots ("all QC-pass samples" and "discovery training samples") | `PLOT_GROUPS`, `:~200` | merge, or say that they are the same samples |
| R7 | No warning that group and batch are confounded in the plot/PC sense: each study holds one sample of one group (SRP219837: 1 normal vs 2 tumour; SRP301216: 2 normal vs 1 tumour), PC1 follows the study, and 2 of 6 leave-one-out refits were impossible. The scanner only checks that every batch has both groups | `scripts/input_scan.py` (batch check), `nextflow/bin/cohort_qc.R`, report | add a cohort-QC line "variance explained by study vs group on PC1/PC2" and a warning when a batch has a single sample of a group |
| R8 | The service's validation of an empty input says only "no FASTQ pairs found"; it does not say that `metadata.csv` or FASTQ files are missing in the folder | `scripts/input_scan.py` message, shown by the service | name the likely causes in the message |

### Accuracy upgrades suggested by this check (to be ranked after the research in section H)
| ID | Upgrade | Type |
|---|---|---|
| A1 | Make rule 3a robust to Salmon variation (median of several Salmon runs, a margin, inferential replicates, or featureCounts-only calls with Salmon as a flag); see U5 | **scientific** |
| A2 | Report batch/group confounding strength and PCA variance by study vs group (R7); refuse or warn below a minimum number of samples per group per batch | technical (warn) / **scientific** (refuse) |
| A3 | Paired design (patient term) when `patient_id` repeats | **scientific** |
| A4 | Add `tests/regression/compare_runs.py` to the release checklist: every new version is compared with a baseline run before it is called validated | technical |
| A5 | Sample-identity QC (sex check against metadata, swap/contamination flags) and richer RNA QC metrics | technical (flags) |

## H. Research findings (2026-10-22; all four research agents have reported)
Evidence levels used below: **V** = a page was opened and read; **S** = seen only in a search-result summary; **M** = recalled from memory, not verified. The agents opened few full texts (some pages failed to load), so every item must be re-read in its original source before it is quoted in a report, a paper or a capability statement. Agent statements are inputs, not conclusions.

### H1. Salmon non-repeatability (agent 1; one source re-opened by me)
- **Why (V, maintainer statements on GitHub issues #613 and #48):** reads are streamed by a background thread, so read order, and with it the start and path of the online/EM optimisation, can change; transcripts with high inferential uncertainty then change. The maintainer calls this intended. One reporter saw variation on a single core, as we did (our test: `-p 1` vs `-p 1` differ by 2.6 % of reads). An older reply (#48) says a seed would only help with one thread; that conflicts with the single-thread observations and with our measurement, and our own test is the stronger evidence.
- **No seed or determinism option in Salmon 1.10.x (V, official docs; also our `--help-verbose`).** `--skipQuant` and `--rangeFactorizationBins` are not fixes.
- **A newer Salmon line claims determinism (V, release page https://github.com/COMBINE-lab/salmon/releases/tag/v2.6.0, opened by me):** "makes deterministic quantification the default"; `quant.sf` "byte-identical at every thread count, on real data"; `--gcBias`, `--seqBias`, `--posBias` supported. Not stated on that page: the release year, the implementation language ("port"), `--validateMappings`, index compatibility with our 1.10.3 index, equivalence of the estimates with 1.10.3. **Must be validated before any use; it would change our numbers.**
- **Maintainer-endorsed route for uncertainty (V):** inferential replicates (`--numGibbsSamples` or `--numBootstraps`, 20-30) with swish (Bioconductor `fishpond`).
- **nf-core/rnaseq (V):** default `star_salmon`; the pseudo-aligner route is optional; no reproducibility advice beyond pinning versions.

| ID | Upgrade (replaces/refines A1, U5) | Changes results | Needs approval |
|---|---|---|---|
| H1a | Treat Salmon as a flagged cross-check; featureCounts stays the primary call (already the primary DE route) | yes (borderline candidates are flagged) | yes |
| H1b | Test Salmon 2.x: run the demo samples twice, require byte-identical `quant.sf`; compare with 1.10.3 within a stated tolerance; check index and flag support | possibly | yes |
| H1c | Loosen rule 3a on the Salmon side to a concordance test (same direction, looser padj or effect size) | yes | yes |
| H1d | Add Gibbs replicates and a swish or uncertainty-aware stability filter on candidates only | yes (modest) | yes |
| H1e | Record Salmon's non-determinism and pin versions in every report | no | no |
| (not recommended) | Averaging a fixed number of Salmon runs: no source, no variance estimate | marginal | - |

### H2. Differential expression with few, partly confounded samples (agent 2)
- **Power (V, Schurch 2016, PMC4878611; yeast, 48 replicates):** at 3 replicates 9 of 11 tools found only 20-40 % of the DE genes (over 85 % for changes above 4-fold); DESeq2 and edgeR kept FDR at or below about 5 % at all replicate numbers; the authors recommend at least 6 replicates, and at least 12 to find most DE genes. Ching 2014 (V, PMC4201821): DESeq2 and edgeR "tend to give the best performance in general", a minimum of five samples in their setting, and paired designs "significantly enhance" power. Our 3-versus-3 demo is underpowered, not just noisy.
- **Pre-filter and shrinkage (V, DESeq2 vignette):** "smallest group size" pre-filter matches ours; apeglm is the recommended shrinkage.
- **Confounding (S/M, Bioconductor support threads, Nygaard 2016):** when group is nested in batch only the group term is estimable and batch cannot be corrected; adjusting for batch when groups are unbalanced across batches can create false group differences. In our demo each study holds both groups, so the design is estimable, but each cell has one or two samples and PC1 follows the study. Do not feed ComBat-seq or SVA/RUV output back into a model that also contains batch (M).
- **Leave-one-out (S):** no published standard procedure and **no source for the 90 % threshold**; related work uses delete-d jackknife or bootstrap stability (AUCOR). With 10 samples one refit moves a gene's pass rate by 10 points, so the rule is effectively "must pass every refit".

| ID | Upgrade | Changes results | Needs approval |
|---|---|---|---|
| H2a | Warn "exploratory, underpowered" below 6 samples per group and when any batch-by-group cell has one sample; report the batch-by-group table, design rank and PCA by study and by group (extends R7, A2) | no | no |
| H2b | Estimability gate: stop when the model matrix is not full rank or group is nested in batch | yes (refuses some analyses) | yes |
| H2c | Paired design `~ patient + group` when `patient_id` repeats; drop unpaired patients in a paired analysis (A3) | yes | yes |
| H2d | Report "N of M refits feasible"; when fewer than 10 are feasible use batch- and group-stratified subsampling or rank concordance instead of skipping | yes | yes |
| H2e | Recalibrate the 90 % leave-one-out threshold and other inherited thresholds on more than one cohort | yes | yes |
| H2f | Cross-check candidates with edgeR or limma-voom | yes | yes |
| H2g | Guard: refuse SVA/RUV/ComBat-seq adjustment when group and batch are confounded | no | no |

### H3. Honest model performance, cutoff and novelty (agent 3; search summaries only, no full texts)
- **Bias of our in-sample numbers (S):** feature selection outside the CV loop makes error estimates strongly optimistic (Ambroise & McLachlan 2002, PNAS); nested CV is the remedy (Varma & Simon 2006, BMC Bioinformatics); bootstrap optimism correction must repeat the whole selection in every resample (Harrell; Steyerberg 2001); stability selection and per-gene selection frequencies address the 8-versus-9-gene problem (Meinshausen & Buhlmann 2010).
- **Cutoff (S):** an in-sample optimal cutpoint inflates performance and is unstable (`cutpointr`, arXiv:2002.09209); report a bootstrap CI for the cutoff.
- **Reporting (S):** TRIPOD+AI (BMJ 2024), STARD 2015, REMARK; MAQC-II as the community reference for predictive-model practice. Regulatory remarks from news summaries (the FDA's laboratory-developed-test rule being vacated in 2025; FDA treats a product as an IVD by intended use) are **not verified**; counsel must confirm.
- **Novelty (no conclusion possible):** no published method matching the fold-median Youden cutoff or the unweighted up/down score was found, but the search was shallow (about four queries). Closest: `cutpointr`, `blkbox` (nested-CV biomarker package, bioRxiv), the sepsis MetaScore (up-minus-down gene score, M), `singscore`. **Do not claim novelty.** Describe the work as documented, pre-specified and reproducible, and run a systematic search (PubMed, Embase, CRAN/Bioconductor, GitHub) first.

| ID | Upgrade | Changes results | Needs approval |
|---|---|---|---|
| H3a | Label in-sample metrics "apparent, optimistic" everywhere (partly done) and add the limitations text | no | no |
| H3b | Fix and log the candidate order and the seeds so the panel is reproducible (does not make it stable) | yes (picks one panel) | yes |
| H3c | Nested or repeated CV with selection and cutoff inside every fold; bootstrap optimism correction repeating the whole pipeline | yes (new honest metrics) | yes |
| H3d | Stability selection: per-gene selection frequencies and a consensus panel | yes (may change the panel) | yes |
| H3e | Bootstrap CI for the cutoff; compare the median-of-folds rule with alternatives | no unless the rule changes | yes if changed |
| H3f | Compare lambda.min and lambda.1se; add the weighted score as a sensitivity analysis | no | no |
| H3g | TRIPOD+AI, STARD and REMARK checklist in the report | no | no |
| H3h | External validation with a pre-specified sample size | yes | yes |
| H3i | Systematic novelty search before any claim | no | yes |

### H4. Landscape, extra QC and licensing (agent 4; plus two licence files I opened myself)
**Licence of the exact versions we run (V, version-tagged LICENSE files opened by me):** STAR 2.7.11b = **MIT** (`alexdobin/STAR` tag 2.7.11b, copyright 2019 Alexander Dobin); Salmon v1.10.3 = **GPLv3** (`COMBINE-lab/salmon` tag v1.10.3). The bioconda package metadata on this pod says GPL-3.0-or-later for both, which is therefore **wrong for STAR** (and it is the source of my earlier statement that STAR is GPL). The agent reported that Salmon's current master `LICENSE` is BSD-3-Clause; that holds for newer code, **not for 1.10.3**. A newer Salmon line (v2.x, the one that claims deterministic output, H1b) may therefore also change the licence situation: check its own tag's `LICENSE` before relying on it. Other licences (package metadata and CRAN/Bioconductor pages): subread/featureCounts GPL-3 (the Rsubread package page says GPL >= 3; standalone subread is M), glmnet GPL-2, pROC GPL >= 3, DESeq2 LGPL >= 3, msigdbr package MIT (the package page states no data terms), fgsea MIT, clusterProfiler Artistic-2.0.

**Hosted use (V, GNU GPL FAQ, first 100k of 156k characters read):** the GPL does not require publishing modifications that only run privately as a web service (the AGPL does); pipes, sockets and command-line arguments usually indicate separate programs. **Distribution (shipping a Docker image, CLI or on-prem install) is different and triggers the GPL terms** (M). Bundling GPL R packages inside distributed R code would make the combined work GPL (M). Counsel must confirm all of this.

**Data and gene-set terms:** MSigDB is CC BY 4.0 for v6.0-v7.5.1 and v2022.1 onward with some sets carrying extra restrictions; the page is silent on versions between; Hallmark and GO are not named on it; KEGG MEDICUS is CC BY-SA 4.0 (V). GO is CC BY 4.0 with attribution (V). **KEGG needs a commercial licence for non-academic use (V); keep it off, including KEGG sets inside msigdbr.** GENCODE: the data-access page states no licence and points to EBI terms (no additional restrictions beyond the data owners' terms, attribution expected; commercial terms not stated) (V). GEO/SRA: NCBI places no restrictions, but submitters may claim IP (V). TCGA/GDC: commercial-use terms not found (open item, as in `docs/PIPELINE_SPEC.md:181`).

**Landscape (V for what the pages state, otherwise not stated):** nf-core/rnaseq runs no significance tests and points to nf-core/differentialabundance, which offers DESeq2 or limma-voom, optional GSEA and a report, with no biomarker or predictive modelling mentioned; Galaxy tutorials cover STAR, featureCounts and DESeq2; Partek Flow is sold by quote (price not public); Qlucore: quote on request; Illumina ICA/DRAGEN lists per-gigabase licence costs. **No product found bundles differential expression, signature discovery and nested-CV validation, but the survey was not exhaustive**, so differentiation is a hypothesis.

**Extra QC the sources describe (thresholds only where a source gives them; none are published pass/fail cut-offs for these metrics):** strandedness check (RSeQC `infer_experiment`; nf-core infers it); gene-body coverage / 3' bias (RSeQC, Picard definition); read duplication, TIN, read distribution (RSeQC); rRNA and intronic fraction (Picard `PCT_RIBOSOMAL_BASES`, `PCT_INTRONIC_BASES`); replicate correlation and depth (ENCODE 2011: replicates below 0.9 correlation should be repeated or explained; about 30 M paired-end reads for polyA+); sex check with XIST and Y genes against metadata, with cut-offs derived from the cohort's own bimodal distribution (informal sources only); swap and contamination (NGSTroubleFinder). RIN is metadata to record, not a computed check.

| ID | Upgrade | Changes results | Needs approval |
|---|---|---|---|
| H4a | Strandedness check against the library type used for featureCounts (the pipeline already infers `-s` from Salmon: verify it independently) | yes if a wrong `-s` is currently used | yes |
| H4b | Sex check against metadata, flag-only at first; sample-swap flag | flag: no; exclusion: yes | flag no / exclusion yes |
| H4c | Add rRNA %, intronic %, 3' bias or TIN, duplication, replicate correlation to the QC table with cohort-derived outlier flags (no fixed pass/fail until calibrated) | flags: no | no (flags) / yes (gates) |
| H4d | Reuse MultiQC / nf-core QC modules for reporting only | no | no |
| H4e | Keep KEGG off; pin and print MSigDB, GENCODE and tool versions and attributions in every report | no | no |
| H4f | Licence and ownership review before any hosted or shipped product: counsel to confirm the points above, add a `LICENSE`, settle ownership of AI-assisted code, check data-use terms of every development cohort | no | yes |
| H4g | Check the v2.x Salmon tag's `LICENSE` and behaviour together with H1b | possibly | yes |

### G2. Numerical reproducibility depends on BLAS threads (found while reviewing the upgrade agents, measured by me, 2026-10-22)
The 19-sample full-cohort DESeq2 refit reproduced the primary table to 5e-14 in yesterday's validated run, but 4.6e-5 (shrunk log2FC) when an upgrade agent ran it with limited threads. Controlled test on the same inputs (`deseq2.R`, 19 samples): with default threads my run equals yesterday's pipeline output **exactly** (difference 0); with `OPENBLAS_NUM_THREADS=1` the shrunk log2FC differs by up to 4.6e-5 (unshrunk 5.7e-10, padj 5.8e-11) and with 4 threads by up to 4.0e-5; 1 vs 4 threads differ by 8.5e-5. R uses OpenBLAS 0.3.34 and the pipeline pins no BLAS thread count, so the numbers depend on the CPU count of the machine (summation order). **Effect on calls in this cohort: none** (0 of 15,111 genes change rule-1 status at 1 or 4 threads; 10 genes lie within 0.001 of the |shrunk log2FC| = 1 cutoff with padj < 0.06, none within 1e-6 of padj = 0.05), but a different cohort or machine could flip a knife-edge gene.
| ID | Upgrade | Changes results | Needs approval |
|---|---|---|---|
| G2a | Pin the BLAS/OpenMP thread count for every R step (for example `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` in the process environment) so results do not depend on the machine; then re-baseline (changes the numbers at the 1e-5 level) and re-run `compare_runs.py` | yes, about 1e-5 | yes |
| G2b | Record the BLAS library and thread setting in the manifest and report | no | no |
| G2c | Until pinned: compare runs from different machines with a tolerance (1e-4 on shrunk log2FC) and list genes within 1e-3 of a cutoff | no | no |

### I. Upgrade batch (2026-10-22): **APPLIED to the live tree on the owner's approval**; backup of the replaced files and the procedure are recorded in `/workspace/Biomarker_pipeline_upgrade/UPGRADE_SUMMARY.md`; new hash list `docs/FROZEN_FILES.upgrade_20261022.sha256`
Four agents implemented the no-result-change items and the opt-in scientific options in `/workspace/Biomarker_pipeline_upgrade/code`; the combined diff, tests and the decisions needed are in `/workspace/Biomarker_pipeline_upgrade/UPGRADE_SUMMARY.md` (diff: `UPGRADE_COMBINED.diff`). Covered ids: R1-R6, U25, R8, H2a (warnings), H2b (estimability gate), H2d (feasibility file), H1a/H1c (opt-in concordance mode, default strict), H3a-H3f (internal validation step), H4b (flag-only sex check), H4e (versions block), H4f (licence table, counsel to confirm), A4 (`compare_runs.py`). Not done: H1b (needs permission to install Salmon 2.x), H2c paired design, H2e threshold recalibration, H2f edgeR cross-check, H3h external validation, H3i novelty search, G2a BLAS pinning (needs approval).

### H0. Ranked plan from sections G and H (my synthesis)
1. **No change of results, do first:** report fixes R1-R8; warnings and tables H2a; labels H3a; `compare_runs.py` in the release checklist (A4); record Salmon non-determinism (H1e); checklist H3g; versions and attributions in the report (H4e); QC flags without gates (H4c, H4d).
2. **Needs your approval, biggest accuracy gain:** the Salmon decision (H1a-H1c, test H1b first); nested CV, stability selection, cutoff CI (H3c-H3e); the estimability gate (H2b).
3. **Larger work:** paired design (H2c), re-calibrated thresholds (H2e), independent validation (H3h), systematic novelty search (H3i), licence and ownership review (H4f) before anything is hosted or shipped.
