# Integration guide: one-command RNA-seq biomarker discovery pipeline

For the team building the user interface. The backend is a single Nextflow command; the UI only has to place files in a folder, start the command and read the output folder.

## 1. Input folder (one folder per job)
`metadata.csv` and the FASTQ files go in the same folder. `metadata.csv` is found automatically as `<input_dir>/metadata.csv`; no other parameter is needed.

### FASTQ files
- Paired-end only, gzip-compressed, extension `.fastq.gz` or `.fq.gz` (case as written; other extensions such as `.fastq` are ignored, and a sample that only has such files is reported as "metadata row without FASTQ").
- One R1 and one R2 per sample. Accepted mate naming (tried in this order, `conf/input.config` `fastq_mate_tags`):

  | Pattern | Example files | Sample ID |
  |---|---|---|
  | `_R1_001` / `_R2_001` (Illumina) | `S1_R1_001.fastq.gz`, `S1_R2_001.fastq.gz` | `S1` |
  | `_R1` / `_R2` | `S1_R1.fastq.gz`, `S1_R2.fastq.gz` | `S1` |
  | `.R1` / `.R2` | `S1.R1.fastq.gz`, `S1.R2.fastq.gz` | `S1` |
  | `_1` / `_2` (SRA) | `S1_1.fastq.gz`, `S1_2.fastq.gz` | `S1` |

- **One naming style per folder.** Mixing styles in one folder (for example some samples `_R1/_R2` and others `_1/_2`) stops the run at the input check and lists the files.
- **One pair per sample: merge lanes first.** Files of several lanes of one sample (`S1_L001_R1_001`, `S1_L002_R1_001`) are not merged: each lane would be a separate sample with the ID `S1_L001`, `S1_L002`, and every ID would need its own metadata row. Concatenate the lane files (all R1 into one R1, all R2 into one R2, same order) into one pair per sample before the run.
- Sample IDs may contain letters, digits, `.`, `_` and `-`, and must start with a letter or digit.
- Input errors are all reported together and stop the run before any computation: R1 without R2 (or the reverse), the same sample and mate given twice (for example `.fastq.gz` and `.fq.gz`), broken links, empty or non-gzip files, unsafe IDs, metadata rows without FASTQ, FASTQ without a metadata row.

### metadata.csv
- Comma-separated, header in the first row, column names exactly as below (case-sensitive).
- **Required columns: `sample_id` and `group`.**
- **Optional columns: `batch_or_study` and `patient_id`** (also optional: `batch`, `study`, `label`, any other columns are carried along).
  - `batch_or_study` missing or empty: all samples are placed in one batch and no batch adjustment is done (a NOTE is printed). If it is given, the design is `~ batch + group` and every batch must contain both groups. A misspelt column name is treated as missing, so the UI must write the exact name.
  - `patient_id` missing or empty: every sample is treated as its own patient (paired tumour/normal samples from one patient are then not recognised as pairs).
- `sample_id` must equal the FASTQ name without the mate tag and extension, be unique, and have exactly one metadata row per FASTQ sample.
- `group` values are matched case-insensitively against `group_aliases` in `conf/input.config`; an unknown label stops the run. Exactly two groups are needed, one of them the baseline `normal`. Current aliases (canonical group on the left):
  - `tumour`: tumour, Tumor, tumor, T, cancer, case, disease, diseased
  - `normal`: normal, Normal, N, adjacent normal, control, healthy, ctrl
  - Deliberately not mapped, because their direction is ambiguous: patient, positive/negative, treated/untreated, wt/mutant. Add such labels in `conf/input.config` only after deciding which side is the baseline. The outputs always show the canonical names (`tumour`/`normal`); the mapping is logged in `00_input/group_mapping.tsv`.

### What cohort size each stage needs (the pipeline accepts any size; stages that cannot be done are skipped with a warning)
Discovery is re-run on every dataset: the differential expression, the candidate rules (including the 19... leave-one-out refits), the LASSO gene selection, the score and the cutoff are all recomputed from the samples in the current `input_dir`. No gene list from an earlier cohort is stored in or read by the discover mode (the locked nine-gene model is used only by `--mode apply`).

| Stage | Needs | When the cohort is too small (or too large) |
|---|---|---|
| Input check | 2 groups (one the baseline), at least 3 samples per group, every batch contains both groups, matching FASTQ and metadata | The run stops before any computation, with the reason. |
| QC, counts, DESeq2, pathway analysis | at least 3 samples per group | As above (these stages always run when the input check passes). |
| Candidate rules (rule 3b, leave-one-out robustness) | at least 10 refits, so at least 10 samples (fewer if a refit would leave a batch with one group) | A WARNING is printed at the start and in the report. The candidate table is still written, but no candidate is declared ("RULE 3B NOT EVALUABLE"). |
| LASSO gene selection, score, cutoff | at least 2 candidates that passed all rules | "NO MODEL FITTED" is written to the model files and the report, the run still finishes with exit code 0. Quantification, differential expression and pathway results are kept. |
| Leave-one-out with many samples | at most 60 samples (`robustness_loo_max_samples` in `config/questions.yaml`) | Above 60 the run fails at the discovery input check, after quantification, so no report is written (`--stop_after analysis` ends the run cleanly before the discovery stage). Switching to the subsample method is a scientific decision for the owner and is not done automatically. |

A small cohort therefore still gets QC, counts, differential expression and pathways, but no gene panel. "No model fitted" never means that no biomarker exists. With the 19-sample cohort all stages ran; whether a panel found in a cohort of this size is stable is a separate question (the nine genes changed with the feature order in the audit).

## 2. Command
```
nextflow run nextflow/main.nf --input_dir <folder> --output_dir <new results folder> --mode discover
```
- `--mode apply --cohort development|validation` scores new samples with the frozen model in `lock/diagnostic_v1` (that lock was built on colon tumour vs normal: use it only for that disease and never for other tissues).
- Optional: `--pod_cpus N --pod_mem_gb N` (default: detected from the container; CPUs = container CPUs, memory = cgroup limit minus 6 GB), `--stop_after quant|analysis|discover`, `-resume`.
- Nothing else needs editing. `nextflow.config` ships with `input_dir`/`output_dir` = null.

### Checking an input before a long run
- `--check_only true` runs only the startup checks: software, reference and index files, memory, disk, and that the input folder has at least one FASTQ file and a `metadata.csv`. It does **not** check pairing, file names, metadata columns, group labels or sample counts, and it still creates the output folder with `logs/`.
- To validate the real input (pairing, naming style, metadata, labels, minimum 3 per group) without computing anything, run the scanner on its own; it writes only to a new, non-existing `--out-dir` and exits 2 with all problems listed if the input is invalid:
  ```
  /workspace/Biomarker_pipeline/env/prefix/pyqc/bin/python3 -I /workspace/Biomarker_pipeline_onewf/scripts/input_scan.py \
    --input-dir <input_dir> --mode discover --out-dir <new empty path> \
    --config /workspace/Biomarker_pipeline_onewf/conf/input.config --batch-optional --patient-optional \
    --min-per-group 3 --questions /workspace/Biomarker_pipeline_onewf/config/questions.yaml
  ```
  The same checks run automatically at the start of every real run, before any computation.

## 3. Output folder
Everything is written under `<output_dir>`:
- `00_input/` (scan, resolved metadata, group mapping), `01_fastp/`, `02_alignment/`, `03_quant_salmon/`, `04_counts/` (count matrix), `05_cohort_qc/` (QC table), `06_dge/`, `07_pathway/`, `08_multiqc/`, `09_discover/{00_input,01_de,02_candidates,03_model}/` (apply mode: `09_locked_model/`), `contract/`.
- **Report: `10_report/report/report.html`** and `report.md`, with `10_report/tables/`, `10_report/plots/`, `10_report/manifest.json`, and the candidate model in `10_report/model/diagnostic_v2_candidate_<run>/`.
- `logs/` (trace, timeline and report per launch), `run_info.json`, `params.json`, `sample_map.tsv`, `tool_versions.tsv`, `RESULTS_INDEX_RUN.md`, `work/`.

## 4. Status for the UI
- Exit code 0 = finished; non-zero = stopped. The reason is in `logs/` and, for input problems, in `00_input/problem_files.tsv`. Do not read the report header's "completed" as proof of completion; use the exit code.
- Quality flags are in `04_counts/sample_qc.tsv` (PASS/WARN/FAIL; the cohort-level QC plots and flags are in `05_cohort_qc/`). FAIL samples are excluded automatically and listed.

## 5. Parallelism and resources
Per-sample steps run in parallel as far as RAM allows (STAR needs about 33 GB and Salmon about 28 GB per job; the number of parallel jobs is derived from the pod memory). From DESeq2 on, all samples run together. Minimum pod: 32 GiB RAM (STAR route); validated on 16 CPUs and 64 GB RAM. The startup check estimates about 13 GB of work space per sample plus 5 GB and warns if the reported free space is lower (the shared volume may report cluster space, not your quota).

## 6. Fixed parameters (do not expose in the UI)
Thresholds are in `conf/thresholds.config`; candidate rules in `config/candidate_rules.yaml`; model seeds, folds in `config/questions.yaml`. They are the values of the validated 19-sample run.

## 7. Known limits (do not hide from users)
- Salmon quantification is not exactly repeatable run to run for some samples (about 2% of reads in paralog groups); the featureCounts route is. It can move a few genes that sit at the rule-1 cutoffs.
- The nine-gene model is validated for reproducibility only; it has no independent validation. Reports must not claim clinical validity.
- Groups with fewer than 10 samples are labelled "MECHANICS TEST ONLY" in the report. A panel found in a small cohort is a hypothesis to be validated, not a result.
