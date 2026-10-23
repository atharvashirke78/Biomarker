# Folder entry: one command from a data folder (2026-10-19, development; not committed)

Plan: `results_rc1prep_19/phase1_plan_20261016/PHASE1_PLAN.md` section 7 and owner decisions 9a. Settings: `conf/input.config`.

## Input folder

```
DATADIR/
  S1_R1.fastq.gz  S1_R2.fastq.gz      (or _R1_001/_R2_001, .R1/.R2, _1/_2; .fastq.gz or .fq.gz; one naming style per folder)
  S2_R1.fastq.gz  S2_R2.fastq.gz
  metadata.csv                        sample_id, group, batch_or_study [, patient_id, label, age, sex, fastq_R1, fastq_R2, ...]
```

- `group`: Tumor / tumor / T / cancer / case / disease / diseased -> tumour; Normal / N / adjacent normal / control / healthy / ctrl -> normal (case-insensitive, whole string; `conf/input.config` `group_aliases`). Every mapping
  used is written to `input/group_mapping.tsv` and the report. Any other label stops the run with the allowed list.
- `batch_or_study` fills both batch and study; separate `batch` / `study` columns override it.
- `patient_id` is needed for discover (patient-level hold-out and input checks); it is never used in a score.
- **This page describes the `scripts/rnaseq.sh` route, where `batch_or_study` is required and `patient_id` is required for discover** (the scanner is called without `--batch-optional` / `--patient-optional`, `scripts/rnaseq.sh:74`). In the one-workflow route (`nextflow run nextflow/main.nf --input_dir ...`, see `docs/README_one_workflow.md`) only `sample_id` and `group` are required: both other columns are optional (`modules/folder_input.nf:29-30`) and its output layout is `OUTPUT/10_report/report/report.html`, not the run folder described below.

## Commands (project root)

```
scripts/rnaseq.sh --make_metadata --input_dir DATADIR                       # draft metadata from the FASTQ names, then stop
scripts/rnaseq.sh --mode quant    --input_dir DATADIR [--outdir OUT]
scripts/rnaseq.sh --mode discover --input_dir DATADIR [--outdir OUT]        # QUANT -> DISCOVER (diagnostic)
scripts/rnaseq.sh --mode apply    --input_dir DATADIR --cohort validation   # QUANT -> locked model lock/diagnostic_v1 (hash-checked)
# counts-level (no FASTQ step): add --counts counts.tsv [--counts_secondary salmon_counts.tsv --rds_secondary tximport.rds] [--qc_table qc.tsv]
# resume an interrupted run:  same command + --run_id <run folder name> --resume
```

## What happens

1. Folder scan (`scripts/input_scan.py`): the folder is listed with lstat. The run stops and lists every problem file: broken link,
   empty file, not gzip, no mate tag, mixed naming styles, R1 without R2 (or R2 without R1), the same sample twice, ID not matching
   `sample_id_regex`.
2. metadata.csv checks: required columns, duplicate IDs, a metadata row without data (stop), data without a row (stop in quant/discover;
   group `unknown` and listed in apply), group labels, empty batch, numeric age, patient_id for discover. No metadata.csv in
   quant/discover: `metadata.DRAFT.csv` is written into the run folder and the run stops. Nothing is written into DATADIR.
3. Run folder `OUT/<mode>_<YYYYMMDD_HHMMSS>/` (OUT defaults to `DATADIR/result`), never overwritten:
   `run_info.json` (commit, uncommitted files, lock check, pod, attempts), `metadata.used.csv`, `input/`, `pipeline/`,
   `nf_<stage>/` (Nextflow logs, traces), links `qc/ align/ counts/ sample_qc/ deg/ candidates/ scores/`,
   `model/diagnostic_v2_candidate_<run id>/` (a candidate, never written into lock/), `tables/`, `plots/`, `tool_versions.tsv`,
   `manifest.json` (sha256 of inputs, config, lock files, report files), `report/report.md` + `report.html`.
4. Resume: `--run_id NAME --resume` re-runs Nextflow with `-resume` in the same folder; refused if the folder is missing, was
   started in another mode or on another folder, or metadata.csv changed. A failed attempt writes `report_partial_<time>/`.

## Known limits

- Research use only; not a clinical or diagnostic device.
- Input: human, paired-end, bulk RNA-seq only. Single-end reads, multi-lane FASTQ (L001/L002) and other species are not supported yet.
- Not built yet: `--gene_id_type` (Entrez/symbol counts), a species guard, pathways inside DISCOVER (Reactome/GSEA port), the cohort
  outlier flag for small APPLY cohorts (single-sample baseline DEFERRED).
- The locked model `lock/diagnostic_v1` is valid only for colon tumour vs normal tissue, scored on featureCounts counts with GENCODE v50 IDs.
- Memory: STAR with the current index (`ref/star_index`, about 26.5 GB loaded) needs a pod with enough real memory. Measured failure:
  STAR was OOM-killed on an 8-CPU / 29.8 GiB pod (2026-10-19 smoke test, `docs/run_log.md`). Recommended: 64 GB RAM and 16 CPUs (the
  configuration of the successful 2026-10-08/09 development run); the strict minimum was not measured.
- The memory guard (`min_mem_star_gb` in `conf/thresholds.config`, read as GiB) stops the QUANT path, also `--stub`, on a pod below that
  limit; it is not bypassed by the pipeline. The FASTQ path of quant/discover/apply has not yet completed on real data (sheet parsing
  tested separately: `tests/integration/parse_sheet_test.nf`).
