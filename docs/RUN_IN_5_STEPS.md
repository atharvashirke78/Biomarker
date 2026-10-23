# Run the biomarker pipeline in 5 steps

Research use only. This pipeline takes paired RNA-seq FASTQ files from colon tissue (tumour vs normal) and gives you gene counts, quality checks, differential expression, pathways, biomarker candidates and a report.

## 1. Put your files in ONE folder
```
my_run/
  S1_R1.fastq.gz   S1_R2.fastq.gz        one pair of files per sample (_R1_001/_R2_001, _R1/_R2, .R1/.R2 or _1/_2; .fastq.gz or .fq.gz)
  S2_R1.fastq.gz   S2_R2.fastq.gz
  ...
  metadata.csv                             one row per sample
```
`metadata.csv` (comma-separated, first row = column names):
```
sample_id,group,batch_or_study
S1,normal,studyA
S2,tumour,studyA
```
`sample_id` = the file name before `_R1`, using only letters, digits, `.`, `_` and `-` (no spaces, commas or quotes); `group` = `normal` or `tumour` (also tumor, T, N, cancer, adjacent normal, case, disease, diseased, control, healthy, ctrl); only `sample_id` and `group` are required: `batch_or_study` and `patient_id` are optional (without `batch_or_study` all samples count as one batch and no batch adjustment is done; without `patient_id` every sample counts as its own patient). Use at least **3 samples per group**, and keep the FASTQ names in one style. Example folder: `examples/make_example_input.sh --out my_example`.

## 2. Start it (one command)
```
/workspace/Biomarker_pipeline/env/bin/nextflow run /workspace/Biomarker_pipeline_onewf/nextflow/main.nf --input_dir /path/to/my_run
```
Optional: `--output_dir /path/to/results` (default: `my_run/Result`), `--mode apply --cohort development` (score the samples with the locked 9-gene model instead of discovering biomarkers), `--stop_after quant` (stop after the counts and QC; see step 4), `--check_only true` (only check that the software, reference, memory and disk are in order and that the input folder has FASTQ files and a metadata.csv, then stop; it does NOT check the file pairing, the names, the metadata contents or the group labels, and it still creates the output folder with a `logs/` folder). Start it from a folder you keep (the resume cache lives there). For long runs: `tmux new -s run1`, then `nohup <the command> > run1.console 2>&1 &`.
To avoid typing the full path: `export PATH=/workspace/Biomarker_pipeline/env/bin:$PATH`.
You can also write `input_dir` (and `output_dir`) once into the USER SETTINGS block at the top of `nextflow/nextflow.config` and start with `nextflow run main.nf` from the `nextflow/` folder; a flag on the command line wins.

## 3. Check progress
`tail -f run1.console` shows the steps and how many samples are done. First lines: **STARTUP CHECKS**, a list of `[ OK ]`, `[WARN]` and `[MISSING]`: fix every `[MISSING]` and start again. Time: about 10 minutes per sample for the alignment (a 3-sample run took about 25 minutes). Files appear in `my_run/Result/` as steps finish.

## 4. Stop and resume
Interrupted (or you used `--stop_after`)? Run the **same command** again, adding `-resume`, from the **same folder**. Finished steps are not repeated. If you changed the input folder, start a new `--output_dir` instead.

## 5. Find the results (`my_run/Result/`)
| Look here | What you find |
|--|--|
| `10_report/report/report.html` | **start here**: the summary with plots and tables |
| `04_counts/`, `contract/` | gene counts per sample, QC table (which samples passed) |
| `06_dge/` | differential expression: volcano, MA, heatmap, tables |
| `07_pathway/` | pathway analysis (over-representation, gene-set enrichment) |
| `09_discover/02_candidates/candidate_table.tsv`, `03_model/` | biomarker candidates and the LASSO model |
| `09_locked_model/scores/` | apply mode: scores and calls of the locked model |
| `00_input/qc_excluded.tsv` | samples removed by the QC, and why |
| `logs/`, `run_info.json`, `params.json` | what was run, with which settings |

## Common messages
| Message | What it means / what to do |
|--|--|
| `NO INPUT FOLDER SET` | no `--input_dir`: add it, or edit USER SETTINGS in `nextflow.config`. The pipeline never guesses a folder. |
| `... does not exist or is not a folder` | typo in the path, or the disk is not mounted |
| `no *.fastq.gz ...`, `metadata.csv not found` | wrong folder, or the metadata file is missing or misnamed |
| `unpaired` / `without a metadata row` / `without FASTQ` | one mate file is missing, or a name differs between the files and `metadata.csv` |
| `does not match sample_id_regex` | a sample name contains a space, comma, quote or another character that is not allowed: rename the FASTQ files and the `sample_id` |
| `group is confounded with batch` | in discover mode, with more than one batch, every batch must contain both groups (here one batch holds only tumour or only normal): add samples or use one batch |
| `unrecognised group label` | the `group` column has a word the pipeline does not know: use one of the labels listed in step 1 (`conf/input.config` `group_aliases`) |
| `too few samples in group` | each group needs at least 3 samples (also after QC removed some) |
| `RUN STOPPED after quantification` | QC removed so many samples that fewer than 3 remain (or fewer than 3 in a group); the message names them and why |
| `STARTUP CHECKS FAILED ... memory` | the pod has too little memory for the aligner: use a 64 GB pod |
| `... conda environments` / `reference` / `STAR index` missing | the install is incomplete or `BIOMARKER_INSTALL_DIR` points to the wrong folder |
| `already holds a run` | that output folder was used: add `-resume` to continue it, or choose a new `--output_dir` |
| `input folder changed since this output_dir was started` | files or metadata changed after the first start: use a new `--output_dir` |
| `needs --cohort` | apply mode needs `--cohort development` or `--cohort validation` |
