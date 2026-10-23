# One workflow: FASTQ folder in, biomarker list and report out

Research use only; not a clinical or diagnostic device. Branch `feat/one-workflow` (not merged, not validated on real data: see "Status").

## Start

`nextflow` is **not on the default PATH** of this install. Use the full path (the wrapper script of the install), or put it on your PATH first:

```
export PATH=/workspace/Biomarker_pipeline/env/bin:$PATH        # optional: then plain `nextflow` works in this shell

/workspace/Biomarker_pipeline/env/bin/nextflow run /workspace/Biomarker_pipeline_onewf/nextflow/main.nf \
    --input_dir <folder> [--output_dir <folder>] [--mode discover|apply] [--cohort development|validation] [-resume]
```

(after merging, replace the `main.nf` path by the one of the installed checkout). Start it from a folder you keep for this project: Nextflow stores its resume cache and `.nextflow.log` in the folder you start from, and `-resume` must be started from that same folder. No `-profile` flag is needed: conda is switched on automatically when `--input_dir` is given.

Only `--input_dir` and `--output_dir` change from run to run. Reference genome, STAR index, annotation, tool environments, thresholds and the
locked model are fixed in the pipeline install (`nextflow/nextflow.config`, `conf/thresholds.config`, `lock/diagnostic_v1/`) and are never
passed on the command line. `input_dir` and `output_dir` can also be written into `nextflow/nextflow.config` (`params.input_dir`,
`params.output_dir`) instead of the command line.

- `--output_dir` default: `<input_dir>/Result`. The Nextflow work folder is `<output_dir>/work`. The folder `Result/` inside the input folder is never scanned as input.
- `--mode discover` (default): fastp, STAR + samtools + featureCounts, Salmon, QC, MultiQC, DESeq2, pathways, candidate rules, LASSO model, biomarker list, report.
- `--mode apply`: the same quantification, then the **locked** 9-gene model (`lock/diagnostic_v1`, hash-checked) scores every sample, then the report. `--cohort` is required (`development` results are labelled optimistic, not validation).
- Install location: `install_dir` (default `/workspace/Biomarker_pipeline`, or the environment variable `BIOMARKER_INSTALL_DIR`) holds the files that are not in git: `ref/` and `env/`. `code_dir` is the checkout that holds `main.nf`.
- Launch from a fixed folder if you want to resume: Nextflow keeps its cache (`.nextflow/`) and `.nextflow.log` in the folder where you start it. Start the resume from that same folder.

Short guide for a first-time user: `docs/RUN_IN_5_STEPS.md`. Startup checks: the first lines of every run list what is installed and what is missing (`[ OK ]`, `[WARN]`, `[MISSING]`); `--check_only true` prints them and stops.
User settings: `input_dir` and `output_dir` can be set once in the USER SETTINGS block at the top of `nextflow/nextflow.config` (command line wins). An empty or invalid `input_dir` stops the run at the start; the pipeline never falls back to a default folder or to test data.

**Stop and continue:** `--stop_after quant` stops after the counts and the sample QC (a QC report is still written), `--stop_after analysis` (discover mode) after DESeq2 and pathways, `--stop_after discover` after the discovery chain; the default runs everything of the chosen mode. In apply mode `--stop_after quant` stops before the locked-model scoring, so the QC can be checked first. Continue with the same command without `--stop_after` and with `-resume`.

## Input folder

```
input_dir/
  S1_R1.fastq.gz  S1_R2.fastq.gz       one pair per sample
  S2_1.fq.gz      S2_2.fq.gz           (the other mate style, but ONE style per folder)
  metadata.csv
```

- **FASTQ names:** `<sample><mate tag><extension>`. Mate tags `_R1_001/_R2_001`, `_R1/_R2`, `.R1/.R2` or `_1/_2` (tried in this order); extensions `.fastq.gz` or `.fq.gz`. The sample ID is the part before the mate tag. Several lanes of one sample are not merged: each lane would be its own sample, so merge the lanes first. (Settings: `conf/input.config`.)
- **metadata.csv columns:**

| Column | Needed | Meaning |
|--|--|--|
| `sample_id` | yes | must equal the sample part of the FASTQ names |
| `group` | yes | `tumour` / `tumor` / `T` / `cancer` / `case` / `disease` / `diseased` or `normal` / `N` / `adjacent normal` / `control` / `healthy` / `ctrl` (case-insensitive, whole string; every mapping used is logged; any other label stops the run; exactly two groups, one of them `normal`) |
| `patient_id` | no | used for the patient-level checks; never used in a score. If missing or empty, every such sample is treated as an independent patient (paired tumour/normal samples are then not recognised as pairs) and the report says so |
| `batch_or_study` | no | or separate `batch` / `study` columns. Without it all samples go into one batch `batch1`, the design is `~ group`, and the report says so |
| `label`, `age`, `sex`, others | no | carried along |

## What stops the run, before any compute

A clear message and exit; nothing is run. Checked by `scripts/input_scan.py` (tested: `tests/integration/test_input_scan.py`, `tests/integration/one_workflow_stub.sh`):

missing `metadata.csv` (a draft is written to the scan folder) · unpaired FASTQ (R1 without R2) · FASTQ without a metadata row · metadata row without FASTQ · missing `sample_id` or `group` column · unknown group label ·
a sample ID with a space, comma, quote or other character outside letters, digits, `. _ -` (`sample_id_regex` in `conf/input.config`) · group confounded with batch (discover mode, more than one batch: every batch must contain both groups) · too few samples in a group (discover mode; limit `analysis_min_samples_per_group` = **3**, `conf/thresholds.config`, marked PROVISIONAL; checked again after the sample QC) · broken link, empty or non-gzip FASTQ · lock files not matching `lock/diagnostic_v1/sha256sums_final.txt` · memory guard (STAR needs a pod with at least `min_mem_star_gb` GiB; never bypassed).

## Samples that fail the sample QC

Quantification runs per sample; the sample QC (STAR unique %, featureCounts assigned %, Salmon mapped %, fastp) marks each sample PASS, WARN or FAIL. FAIL samples are **excluded** from the cohort steps and listed with the reason in `00_input/qc_excluded.tsv` (and in the report). If fewer than 3 samples remain, or (discover mode) a group has fewer than 3 samples left, the run **stops after quantification** with a message that names the excluded samples and why. The cohort QC needs at least 3 samples that pass the QC.

## Output (everything under `output_dir`)

| Folder / file | Content |
|--|--|
| `00_input/` | folder scan, resolved metadata, group mapping, samplesheet, clinical tables, `input_summary.json` (notes, e.g. batch default) |
| `contract/`, `sample_map.tsv` | counts matrix, QC table, sample sheet handed from quantification to analysis |
| `01_fastp/ 02_alignment/ 03_quant_salmon/ 04_counts/ 05_cohort_qc/ 08_multiqc/` | one folder per stage |
| `06_dge/ 07_pathway/` | DESeq2 on all QC-pass samples, ORA and GSEA (discover) |
| `09_discover/` | candidate rules, **candidate_table.tsv** (biomarker candidates), LASSO model (`03_model/`) |
| `09_locked_model/scores/` | apply mode: scores, calls, evaluation when every sample has a label |
| `10_report/` | **report/report.html**, `report.md` (input, stages, discovery result, plots, methods and versions, limitations and next validation steps), `tables/`, `plots/`, `manifest.json` (sha256 of inputs, config, lock, outputs) |
| `logs/`, `run_info.json`, `params.json`, `tool_versions.tsv`, `RESULTS_INDEX_RUN.md` | Nextflow trace/report/timeline per launch, run record (attempts, lock check at start and end), all parameters, tool versions, folder guide |

## Resume

Same command plus `-resume`, same `output_dir`, started from the same folder. The input scan is repeated and must equal the stored `00_input`; if the input folder changed, the run stops (use a new `output_dir`). A second start without `-resume` on an existing `output_dir` also stops. If the reference had to be built during the first run (not the case in the standard install, where `ref/` exists), a resume re-runs the sample tasks because the reference paths change. Cheap tasks that embed the run name (`MANIFEST_QUANT`, `PACKAGE_MODEL`, `REPORT`) run again on every launch.

## Status and limits

- Tested with `-stub-run` (structure test, no real data): `tests/integration/one_workflow_stub.sh`. The regression suite is run through `tests/regression/run_regression_branch.sh` and judged by the acceptance rule in `tests/regression/check_numeric_tolerance.py`: all decisions identical (`check_decisions.py`) and every numeric value within a relative 1e-5 of a control run of the frozen code on the same pod. Not yet run on real FASTQ through this entry.
- Existing entries (`-entry QUANT|ANALYSIS|DISCOVER`, no entry without `--input_dir`) keep their behaviour; the discover stage folders are only moved under `09_discover/` when `--input_dir` is given.
- Discover mode runs DESeq2 twice on purpose: once on all QC-pass samples (`06_dge`, pathways) and once inside the discovery chain on the training samples (`09_discover/01_de`). No hold-out is used in this entry.
- Apply mode on FASTQ needs a metadata row and a group label for every sample (QUANT needs the baseline group).
- Only the STAR + featureCounts route; the Salmon-only fallback is not supported in this entry. Human, paired-end, bulk RNA-seq; diagnostic question only (prognostic and predictive are not built).
