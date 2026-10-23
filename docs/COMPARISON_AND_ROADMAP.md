# Comparison with public pipelines, and what we can deliver next

Date 2026-10-23. Four research agents read public pages; this file keeps only claims that were opened (V) or that we re-checked against our own code. Items marked **U** were not verified and must be read at the source before they are quoted. Research use only; nothing here is legal or regulatory advice.

## 1. Our inputs and outputs (checked against the code)

| | |
|---|---|
| In | Paired-end human bulk RNA-seq FASTQ (`.fastq.gz` / `.fq.gz`, four mate-tag styles) and `metadata.csv` (`sample_id`, `group` required; `batch_or_study`, `patient_id` optional). Options: `--mode discover\|apply`, `--cohort`, `--stop_after`, `-resume` |
| Out, always | Counts (featureCounts and Salmon), per-sample QC verdict PASS/WARN/FAIL with reason, cohort QC (PCA, correlation, library size, MAD flags, flag-only sex check), MultiQC, HTML and MD report, `run_info.json`, `params.json`, `tool_versions.tsv`, sha256 `manifest.json` |
| Out, discover | DESeq2 tables and plots, ORA and GSEA, candidate table with rule funnel (rules 1, 2, 3a, 3b), LASSO panel, cutoff, development performance (labelled optimistic), internal-validation files (nested CV, bootstrap optimism, cutoff interval, selection frequency) |
| Out, apply | Locked 9-gene scores and calls per sample; evaluation when labels are complete |
| Service | Upload, validate, queue, resume, cancel, delete, live chart, JSON summary, figures, report download |

## 2. Comparison

| Area | nf-core/rnaseq 3.27.0 and nf-core/differentialabundance 2.0.0 (V) | Us |
|---|---|---|
| Samplesheet | `sample, fastq_1, fastq_2, strandedness` (auto allowed) | Folder scan + `metadata.csv`; strandedness inferred from Salmon only |
| QC | RSeQC (7 modules), dupRadar, Qualimap, Preseq, biotype, Kraken2/Sylph contamination, MultiQC flags strandedness mismatch | Mapping-rate and fastp metrics, PASS/WARN/FAIL verdict, MultiQC, cohort PCA/MAD, sex flag. **No independent strandedness, rRNA, 3' bias, duplication or contamination check** |
| Reproducibility | Container profiles (docker, singularity, apptainer, wave), `software_versions.yml`, `params.json`, release pinning | Conda envs, `tool_versions.tsv`, `params.json`, git commit and uncommitted-file list, sha256 of inputs, configs, lock and outputs. **No containers** |
| Differential analysis | Contrasts file with blocking, DESeq2 / limma-voom / DREAM, gProfiler2, decoupler, Shiny app | One contrast, `~ batch + group`, apeglm. No second DE method |
| Verdict and model | None | Per-sample verdict, candidate rules with leave-one-out, LASSO score, locked model, pre-registered rule |
| Tool-specific overlap | | Nothing public found that combines FASTQ to locked-model scoring (**U**: a shallow search; no systematic review) |

Biomarker-modelling tools (V on CRAN pages): `nestedcv` 0.9.0 (MIT, nested CV with filters inside the outer loop), `cutpointr` 1.2.1 (GPL-3, bootstrapped cutpoints), `stabs` 0.7-1 (GPL-2, stability selection), `mlr3` 1.8.0 (LGPL-3). `blkbox` was removed from CRAN (2020) and `MetaIntegrator` archived (2024): do not depend on them. BioDiscML, MAQC-II and SEQC2 details were not verified.

**Corrections to the agents' "we lack" lists** (we already have these): MultiQC, `params.json`, `tool_versions.tsv`, resume in the service, internal nested CV / bootstrap / cutoff interval / selection frequency (upgrade H3), flag-only sex check (H4b).

**Validated against our code:** the biggest honesty gap stays as stated in `internal_validation.R`: the candidate screen (rules 1-3b) is selected on all samples and is not repeated inside the resamples, so every development estimate is still optimistic. Only the locked external test is a claim.

**Arithmetic re-checked:** 0.9^29 = 0.047 and 0.95^59 = 0.049, so 0 errors in 29 (59) samples per class is the smallest design whose one-sided 95% bound excludes 90% (95%). 7 of 8 correct gives a Clopper-Pearson 95% interval of about 47% to 99.7%. The current rule (at most 1 error per class on about 7-8 per class) is a feasibility gate, not evidence of accuracy.

## 3. What we can deliver more (ranked, value against effort)

**A. Low effort, no change of results**
1. **DONE 2026-10-23:** strandedness check independent of Salmon (RSeQC 5.0.6 `infer_experiment`), flag only, in `05_cohort_qc/strandedness_check.tsv` and the report. Real test on the 6 demo samples: all agree with Salmon; counts byte-identical to the earlier run. Cut-offs are nf-core's values applied to the share among determinable reads (our choice; nf-core's exact formula not checked).
2. Export bundle in the service: counts, DE tables, candidate list, model coefficients, scores as CSV in one download.
3. Inline sample-sheet validation messages in the UI (the scanner already produces them).
4. Report additions from TRIPOD+AI / REMARK / STARD (V for TRIPOD+AI and REMARK): intended use, study flow diagram, reference standard and label check, full model specification and freeze hash, exact confidence intervals (Clopper-Pearson) on every performance number, limitations, code availability.
5. Pin BLAS threads and record BLAS in the manifest (upgrade G2).
6. Fix the open upload race and add retention and quota rules in the service.

**B. Medium effort**
7. QC metrics: rRNA %, duplication, 3' bias or TIN, gene-body coverage, contamination screen; cohort-derived flags first, gates only after calibration.
8. Containers (docker / apptainer) and a pinned release tag, so a run is independent of the pod's conda state.
9. Permutation-label null for the cross-validated AUC (show the pipeline does not give high AUC on random labels).
10. Candidate screen inside the resampling loop (a fully nested estimate). Needs a decision: `nestedcv` with a custom DESeq2 filter, or our own wrapper.
11. Service: user accounts with roles, append-only audit log, per-run sharing. Needed before any compliance wording.
12. Contrasts file with blocking and a second DE method (limma-voom) as a cross-check.

**C. Needs decisions or data**
13. Independent validation: at least 30 tumour and 30 normal from at least 3 cohorts not used in development, model frozen first, metrics and CIs pre-specified. The existing 15-sample set can only serve as a feasibility check.
14. Paired design `~ patient + group`, more than two groups, other species, single-end input.
15. Any clinical claim (LDT or IVD): intended use, analytical validation, locked software with change control, quality system. The regulatory landscape is **U** (the 2025 FDA LDT rule was reported vacated; current status not verified): ask regulatory counsel.

## 4. Systematic review results (2026-10-23; four protocol-driven reviews, files in `docs/review/`)

Each review wrote its protocol first, logged search strings and dates, and marked facts V (opened) or S (snippet only). Counts are hand tallies and approximate. Fetched pages were read through a summariser, so V is not word-for-word checking. Not exhaustive: relevance-ranked, partial bioRxiv coverage, closed-source tools not searched.

| Review | Screened / included | File |
|---|---|---|
| 1 RNA-seq pipelines | 202 / 21 | `01_rnaseq_pipelines.md` |
| 2 Signature methods | about 215 / 41 | `02_biomarker_methods.md` |
| 3 QC and reproducibility | about 184 / 33 | `03_qc_reproducibility.md` |
| 4 Reporting and validation | about 40+ / full guideline texts | `04_reporting_validation.md` |

**Findings, re-checked against our repo where marked (R)**
- Mature public alternatives are nf-core/rnaseq 3.27.0 with differentialabundance 2.0.0 and snakemake-workflows/rna-seq-star-deseq2 3.1.2, all MIT. None of the 21 included workflows has a per-sample gate that excludes samples, a locked-model apply mode or output checksums in its docs (S for the absence; not proof).
- Our gaps: no LICENSE file (R, confirmed), no tagged release or DOI, no containers, no independent strandedness or RNA-QC metrics, no paired design (R: `deseq2.R` has no patient term), human paired-end only.
- Only the LASSO, score and cutoff stage is repeated inside the internal resampling (R: stated in `internal_validation.R`). DESeq2 and rules 1-3b are outside, so internal estimates remain optimistic. Selection must sit inside every resample (Ambroise and McLachlan 2002, Varma and Simon 2006, MAQC-II, nestedcv paper: all V).
- The pre-registered rule (at most 1 error per class, 7 normal and 8 tumour) passes 29% of the time at true accuracy 0.80, 69% at 0.90, 90% at 0.95 (re-computed by me). A perfect result on 15 samples still has exact lower bounds of about 0.6.
- **Corrections to earlier text:** (a) Salmon is GPL-3 up to 1.10.x and BSD-3 from 2.0.0 (licence file at each tag, V); we run 1.10.3 (R). (b) SEQC2 concerns DNA sequencing; the RNA-seq reference is SEQC/MAQC-III 2014, and MAQC-II covers microarray models. (c) The LDT rule was vacated on 2025-03-31 and FDA's page shows a restoring rule dated 2025-09-19; nothing later was confirmed, ask counsel.
- Salmon output is not bit-reproducible (maintainer statement, GitHub issue 185, one source): hash inputs, references, versions and parameters; compare `quant.sf` with a tolerance.
- No source states pass/fail thresholds for rRNA %, 3' bias, TIN, duplication or sex calls. Use cohort-derived flags and label them our own choice. Sources do state: ENCODE replicate Spearman above 0.9 (isogenic) and 0.8 (anisogenic); nf-core drops samples below 5% uniquely mapped reads and uses 0.8 / 0.1 strandedness cut-offs.
- Keep the sex check flag-only (a peer-reviewed paper argues sex-check QC discriminates against transgender, gender diverse and intersex people).
- Use STARD 2015 as the primary guideline for the validation and TRIPOD+AI for the model specification; REMARK fits only partly.

**Sample sizes (95% confidence, per class; computed by the review agent, formula checked against an epiR worked example)**

| Design | Per class |
|---|---|
| Precision plus or minus 0.10 around 90% | 35 |
| Precision plus or minus 0.05 around 95% | 73 |
| Exact test 95% against null 80%, pass if at most 2 errors | 30 |
| Exact test 99% against null 90%, pass if at most 1 error | 46 |
| AUC 0.95, half-width 0.05 | 50 |

Recommended target: at least 50 per class from at least 3 independent studies or sites, model frozen first.

**Revised ranked plan**
1. Licence decision, git history, tagged release, DOI (owner decision on the licence).
2. RSeQC `infer_experiment` as a second strandedness estimate, WARN-only; then Picard RNA metrics, RSeQC TIN and gene-body coverage, dupRadar as cohort-flagged metrics. Needs a new conda env (approval).
3. Fully nested internal estimate: DESeq2 and rules inside each outer fold, a label-permutation null, whole-procedure bootstrap. Report it as the headline internal estimate.
4. Optional container profile and container digests in the manifest; record thread settings; tolerance-based comparison for Salmon.
5. Paired design and a contrasts file (scientific decision).
6. Independent validation as in the table above, reported by cohort with exact CIs and a rule whose operating characteristics are stated.
7. Benchmark our counts and DE on the 19 development samples against nf-core/rnaseq.
8. Service: accounts, roles, audit log.

**Not verified by any review:** STAR and OpenBLAS run-to-run determinism, FDA premarket pathways and IVDR classification for this test type, Qualimap and Kraken2 licence variants, whether pminternal, rms or caret can re-run an upstream DESeq2 step per resample, any literature after about September 2026.

## 5. Sources opened (V)
nf-core/rnaseq usage and output docs, release pages (https://nf-co.re/rnaseq/latest/docs/usage); nf-core/differentialabundance (https://nf-co.re/differentialabundance/latest/docs/usage); CRAN pages for nestedcv, cutpointr, stabs, mlr3, blkbox, MetaIntegrator; TRIPOD+AI (https://pmc.ncbi.nlm.nih.gov/articles/PMC11019967/); REMARK (https://pmc.ncbi.nlm.nih.gov/articles/PMC3362748/); MAQC-II (https://pmc.ncbi.nlm.nih.gov/articles/PMC3315840); ENCODE4 bulk RNA-seq standards; Conesa 2016 (https://pmc.ncbi.nlm.nih.gov/articles/PMC4728800/); Seqera pricing and roles pages; Basepair pricing; Terra. Seen only as snippets (S): STARD item text, Buderer 1996 formula, Illumina, DNAnexus, Latch, Partek, OmicsBox, iDEP, Schurch 2016, LDT court ruling.
