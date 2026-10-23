# Freeze record (DRAFT, not committed)

Status: DRAFT 2026-10-21. The final cycle ran on 2026-10-21; after it, the report and stub-test fixes below were made and the stub suite and the apply example were re-run. Stub suite: 88 of 88 PASS. Nothing in this file has been committed.

## 1. Code state and diff summary

- Branch: `feat/one-workflow`
- Base commit: `46c4860169e5601b6c5f593243d1b05dceba55b4` ("Add HANDOVER_diagnostic.md ...")
- The earlier code freeze commit is `495b4e6` (branch `poc-platform`), recorded in `docs/run_log.md`.
- Uncommitted changes against the base: 28 tracked files changed (503 insertions, 107 deletions), 3 renames into `archive/`, and the new files below.

Changed files:
`.gitignore`, `conf/input.config` (one line: `sample_id_regex`),
`modules/{cohort_qc,contract,deseq2,discover,fastp,fastqc,featurecounts,merge_counts,multiqc,pathway,reference,salmon,sample_qc,samtools,star}.nf`,
`nextflow/main.nf`, `nextflow/nextflow.config`,
`nextflow/bin/{candidate_rules.R,fit_diagnostic_model.R,package_model.py}`,
`scripts/{apply_locked_model.R,input_scan.py,make_report.py}`.

Moved to `archive/`: `modules/verify_fastq.nf.bak_pre_identical_check`, `modules/verify_fastq_with_identical_check.nf.pending`, `nextflow/nextflow.config.bak_pre_salmon2`.

New files:
`README.md`, `docs/{README_one_workflow.md,RUN_IN_5_STEPS.md,VALIDATION_LEDGER_DRAFT.md,FREEZE_RECORD_DRAFT.md}`, `docs/img/`, `examples/`,
`modules/{folder_input,lock_score,report}.nf`, `nextflow/bin/{check_qc_min.py,filter_clinical.py}`,
`tests/integration/{folder_input_test.nf,one_workflow_stub.sh,test_loo_refits.sh}`,
`tests/regression/{check_numeric_tolerance.py,run_diagnostic_branch.sh,run_regression_branch.sh}`.

Not touched: `lock/`, `config/`, the validation samples.

Last change of this cycle (model review fix 3, made opt-in):
- New switch `--patient_level_cv` (default `false`) in `nextflow/nextflow.config`, passed by `modules/discover.nf` to `candidate_rules.R` and `fit_diagnostic_model.R`.
- Default = the locked sample-level method: leave-one-SAMPLE-out, 19 refits and `loo_min_pass` 18 for the 19-sample design, folds by sample.
- With the switch on and a repeated `patient_id`: leave-one-PATIENT-out (11 refits, pass >= 10 for the 19-sample design) and folds that keep a patient together.
- The guard "derived loo_min_pass != locked 18" now applies only to the locked (sample-level) path.
- With repeated patient IDs and the switch off, the report says: "Repeated patient_id found; folds and refits are sample-level (the locked method). Use --patient_level_cv for patient-level." (`scripts/make_report.py`; `candidate_rules_method_note.txt` is unchanged so the 19-sample outputs stay byte-identical).
- `candidate_rules.R --plan-only` is a test hook: it also needs the environment variable `BIOMARKER_TEST_PLAN_ONLY=1`, is passed by no workflow code, and is named in no report or method note.

Report fix after the first cycle (`scripts/make_report.py`, `tests/integration/one_workflow_stub.sh`):
- For `--cohort validation` the report no longer carries the development-performance line. It says instead: "These samples were not used to build the locked model. The locked model (lock/diagnostic_v1) was applied without any change. Quality-control exclusions are listed below."
- For `--cohort development` and for discover mode the line "Gene selection was done on all development samples; development performance is optimistic; only independent validation counts." stays.
- A pass-rule line for the validation report is NOT added yet; it comes later with the validation protocol.
- Stub tests: the validation test checks the new text and that "Gene selection was done" is absent; the development test checks the old line is present and the validation text absent; the two fix-3 fold tests now use `--patient-level-cv true`, plus a new default-off case (sample-level folds, no "kept together" text).

## 2. Acceptance rule (owner decision 2026-10-21)

1. Every decision is identical to the locked results: `tests/regression/check_decisions.py`, all checks PASS.
2. Every numeric value of every compared table is within a relative difference of 1e-5 of the CONTROL run (frozen code, same pod): `|a - b| / max(|a|, |b|) <= 1e-5`, with identical column names, row counts, gene/sample order and identical text cells. Script: `tests/regression/check_numeric_tolerance.py --control <control run> --run <run> --out <tsv>`. Control: `results_rc1prep_19/regression_control_20261021`.
3. Byte level against the control for the 19-sample discover outputs (all files except run-specific provenance, version, manifest and summary files): `tests/regression/check_pod_baseline.py` plus a direct `cmp` of every file (`discover_byte_compare.tsv`).
4. Default behaviour of the refit plan: `tests/integration/test_loo_refits.sh` (19 refits, loo_min_pass 18; switch on gives 11 patient-level refits).

Runner: `tests/regression/run_regression_branch.sh <new dir>` (A discover 19 samples, B diagnostic with locked model, C stub default).

## 3. Results of the final cycle

Three runs, in this order. Folders are under `results_rc1prep_19/`. Pod: 16 CPUs, 64 GB cgroup.

1. **Final cycle, 2026-10-21 12:18-13:05 UTC** (`final_cycle_20261021/`, `STATUS.txt`): stub suite, refit test, regression A/B/C with checks, apply and discover examples.
2. **Re-run after the report fix, 14:02-14:30 UTC** (`final_cycle2_20261021/`): stub suite and apply example only.
3. **Stub suite after the final validation wording** (`final_cycle3_stub_20261021*`): stub suite only.

| Item | Result | Run |
|---|---|---|
| Stub suite (`tests/integration/one_workflow_stub.sh`) | **88 passed, 0 failed** (earlier runs: 83 of 86, then 87 of 88; the failures were the opt-in fold tests and the validation-report line, both fixed) | 3 (final) |
| Refit test (`tests/integration/test_loo_refits.sh`) | 16 of 16 PASS. Default: 19 refits, none skipped, loo_min_pass 18. Switch on: 11 patient-level refits, pass >= 10. | 1 |
| Regression A (discover, 19 samples) | exit 0. 19 refits, all run, `loo_min_pass` 18, 1390 candidates. All 32 discover output files (apart from run-specific provenance/manifest/summary files) are byte-identical to the control (`discover_byte_compare.tsv`). | 1 |
| Regression B (diagnostic, locked model) | 48 of 48 checks PASS. The script exits 1 in this run and in the control, because its last step prints "deseq2_results_all.tsv DIFFERS from results_rc1prep_19/deseq2_results_all.tsv" (pre-existing; not caused by this work). | 1 |
| Regression C (stub default) | exit 0 | 1 |
| Decisions (`check_decisions.py`) | 30 of 30 PASS | 1 |
| Numeric rule (`check_numeric_tolerance.py` vs control) | ACCEPTED: 45 IDENTICAL, 3 NEW_OK, 0 out of tolerance; largest relative difference 0 (tolerance 1e-5) | 1 |
| `check_pod_baseline.py` vs the old pod baseline | rc 1, as for the control itself: 24 DIFFERENT for the control, 26 for this run. The differences are DESeq2/GSEA numbers from the older pod (the reason the control and the 1e-5 rule exist) plus text files. The 2 extra files are `DISCOVER_SUMMARY.md` and `RESULTS_INDEX_DISCOVER.md` (wording from the model-review fixes, not numbers). | 1 |
| Apply example (6 samples, `--cohort development`) | exit 0, 64 tasks completed, all 6 calls match the labels. Peak anon memory 30.53 GiB, oom_kill 0. In run 2 (after the report fix) the 5 files in `09_locked_model` (scores, evaluation, evaluation_sens_spec, excluded_samples, run_info) are byte-identical to run 1; peak anon 30.53 GiB, oom_kill 0. | 1 and 2 |
| Discover example (6 samples) | exit 0, 74 tasks completed. Rule 3b NOT EVALUABLE (4 of 6 refits could run, at least 10 needed), 0 candidates, no model fitted, and the report says so. Peak anon memory 30.58 GiB, oom_kill 0. | 1 |

What was NOT repeated after the report fix: the refit test, regression A/B/C and the discover example (run 1 only). The fix changes report text only (`make_report.py`, and the wording of the validation line after run 2 was changed once more, covered by the stub suite of run 3); the apply example was not re-run for that last wording change.

Memory note: in both example runs the sampler's `memory.current` reached about 59.6 GiB, but almost all of that is file cache (reclaimable); anon was about 30.5 GiB and no process was killed.

Note on the comparison checks of run 1: the first attempt (inside the cycle driver) crashed on a relative path because the check scripts must be started from the install root. `check_decisions.py` and `check_pod_baseline.py` were re-run from there (`*_rerun.*` files in `regression/`); the cycle itself was not repeated.

## 4. Known limits

- Gene selection (candidate rules, LASSO) uses all development samples. The development performance is therefore OPTIMISTIC and is not a validation result.
- Patient-level folds and refits are opt-in (`--patient_level_cv`). The default is the locked sample-level method. In the paired 19-sample design the samples of one patient are not independent, and the default method does not account for that; the W2 input warning states it.
- Single-route input (no second quantification route): rule 3a is skipped and rule 2 is modified, so the candidates are not equivalent to the 19-sample method.
- Small cohorts: with fewer than 10 refits that can run (`MIN_REFITS`, provisional), rule 3b is NOT EVALUABLE; no candidate is declared and no model is fitted. A cohort with fewer than the small-cohort limit per group is a mechanics test only (the report says so). The 6-sample examples are in this category.
- The full 19-sample discover run from FASTQ to report has NOT been done in this form. The 19-sample regression starts from the existing count matrices, and the examples use 6 samples with 100,000 read pairs each.
- Several settings are PROVISIONAL (no source), for example `robustness_subsamples`, `robustness_subsample_fraction`, `robustness_loo_max_samples` and `MIN_REFITS`.

## 5. Lock files

`sha256sum -c lock/diagnostic_v1/sha256sums_final.txt`: 7 of 7 OK (checked 2026-10-21 in this checkout, after all changes of this cycle).
sha256 of `lock/diagnostic_v1/sha256sums_final.txt`: `37657d532d47c63cdf7b5d21c3dbc55133a79fb3e752b004eb05d800e3e53a6c`.
Files covered: `lock/diagnostic_v1/{final_model_genes.tsv,final_model_spec.txt,lock_statement.md,sha256sums.txt,step12a_v2_summary.md}`, `config/candidate_rules.yaml`, `config/reference_genes.yaml`.

## 6. Final pipeline script and the one command

Final pipeline script: `/workspace/Biomarker_pipeline_onewf/nextflow/main.nf`

One command a new person runs (details: `docs/RUN_IN_5_STEPS.md`):

```
/workspace/Biomarker_pipeline/env/bin/nextflow run /workspace/Biomarker_pipeline_onewf/nextflow/main.nf --input_dir /path/to/my_run
```

`my_run` holds the FASTQ files and `metadata.csv`. Results go to `my_run/Result` (change with `--output_dir`). Add `--mode apply --cohort development` to score samples with the locked 9-gene model instead of discovering biomarkers.

## 7. Frozen code identity (updated 2026-10-23, after the RNA-QC and nested-validation additions)

Repository `https://github.com/atharvashirke78/Biomarker`, branch `main`. Code commit `c66b1d6`; earlier points: tag `v0.9.0-rc1` = `c5e3ecf` (before the LICENSE file, RNA-QC metrics and nested validation). This repository starts at `a2c8890` and does not contain the history of `/workspace/Biomarker_pipeline` (base commit `46c4860`, lock tag `lock-diagnostic-v1`).

`docs/FROZEN_FILES.sha256` lists the sha256 of the 196 tracked files in `nextflow/ modules/ scripts/ conf/ config/ lock/ tests/ service/ env/` as committed in `c66b1d6`. Check from the repository root: `sha256sum -c docs/FROZEN_FILES.sha256` (0 failures when written, 2026-10-23). sha256 of that list: `d13cef63b87b46569ec7a2f8d5a37218c825a769c72f1dddebaf8a79a1704dee`.
Lock: `sha256sum -c lock/diagnostic_v1/sha256sums_final.txt` gives 7 of 7 OK. The three older manifests are in `archive/freeze_history/` and no longer match the tree.

Tests at `c66b1d6`: stub suite 91 passed, 0 failed, 0 skipped (the byte-identity guard runs from `tests/fixtures/fit_diagnostic_model.46c4860.R`); service tests 25 of 25; real-browser check of the page (`service/tests/browser_check.py`, headless Firefox) 16 of 16; strandedness and RNA-QC unit tests 14 of 14; nested-validation test (`tests/regression/test_nested_full_validation.sh`, small settings, 19-sample matrix plus a synthetic noise case) passed in its author's run; real 6-sample demo runs: 83 of 83 tasks, featureCounts counts byte-identical to the earlier run, sample QC statuses identical (Salmon mapped % differs in the last digit: known Salmon non-determinism). Not run in this repository: the 19-sample regression and numeric acceptance (they need the 19-sample result folders), the nested validation at default settings, and the 15-sample independent validation.

Sections 1-6 above describe the 2026-10-21 cycle and the old checkout; where they name `/workspace/Biomarker_pipeline_onewf` or "uncommitted", this section replaces them.
