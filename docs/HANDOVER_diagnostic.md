# Handover: diagnostic RNA-seq biomarker pipeline (colon tumour vs normal)

Status date 2026-10-19. Repository `/workspace/Biomarker_pipeline`, branch `poc-platform`. Research use only; not a clinical device.
Every number is marked **VERIFIED** (re-read from the named file or command output on 2026-10-19) or **NOT VERIFIED** (estimate or not checked).
Details: `docs/run_log.md`, `docs/README_folder_entry.md`, `results_rc1prep_19/phase1_plan_20261016/PHASE1_PLAN.md`.

## 1. What the pipeline does

You put paired-end RNA-seq FASTQ files and a `metadata.csv` in one folder and run one command. The pipeline checks the folder and the
metadata, then trims reads (fastp), quantifies them two ways (Salmon, and STAR alignment plus featureCounts), applies per-sample and cohort
QC, and writes a counts matrix. In **discover** mode it then runs differential expression (DESeq2), the candidate-gene rules and a LASSO
model fit, and packages a *candidate* model. In **apply** mode it scores each sample with the **locked** 9-gene diagnostic model and, when
labels are present, reports the pre-registered pass rule. Every run gets its own dated folder with a report, tables, plots and a manifest
of file hashes, and an interrupted run can be resumed.

## 2. How to run it

Input folder:
```
DATADIR/
  S1_R1.fastq.gz  S1_R2.fastq.gz      (or _R1_001/_R2_001, .R1/.R2, _1/_2; .fastq.gz or .fq.gz; one naming style per folder)
  metadata.csv                        sample_id, group, batch_or_study [, patient_id, label, age, sex, ...]
```
- `group` labels such as Tumor, tumor, T, cancer, case, disease, diseased -> tumour; Normal, N, adjacent normal, control, healthy, ctrl -> normal (case-insensitive). Every mapping used is
  logged; any other label stops the run with the allowed list.
- `batch_or_study` fills both batch and study unless separate `batch` / `study` columns exist.
- `patient_id` is needed for discover mode when started with `scripts/rnaseq.sh`; it is never used in a score. In the one-workflow route (`nextflow run nextflow/main.nf --input_dir ...`) only `sample_id` and `group` are required; `batch_or_study` and `patient_id` are optional (`modules/folder_input.nf:29-30`).

One command (project root):
```
scripts/rnaseq.sh --mode quant    --input_dir DATADIR [--outdir OUT]
scripts/rnaseq.sh --mode discover --input_dir DATADIR [--outdir OUT]
scripts/rnaseq.sh --mode apply    --input_dir DATADIR --cohort validation [--outdir OUT]
scripts/rnaseq.sh --make_metadata --input_dir DATADIR        # writes a metadata draft from the FASTQ names, then stops
# resume an interrupted run: same command + --run_id <run folder name> --resume
# counts-level runs (no FASTQ step): add --counts counts.tsv
```
Output of `scripts/rnaseq.sh`: `OUT/<mode>_<YYYYMMDD_HHMMSS>/` with `report/report.html`, `tables/`, `plots/`, `manifest.json`; never overwritten. (The one-workflow route writes `OUTPUT/10_report/report/report.html`, `10_report/tables/`, `10_report/plots/`, `10_report/manifest.json`.)
All cutoffs and settings come from `conf/thresholds.config`, `conf/input.config` and `nextflow/nextflow.config`.

## 3. What is locked

| Item | Value | Status |
|--|--|--|
| Lock folder | `lock/diagnostic_v1/` | VERIFIED |
| Git tag | `lock-diagnostic-v1` (tag object `7d224a8`) -> commit `10bd636`, 2026-10-14 12:32 UTC | VERIFIED |
| Model | balanced up/down score: mean log2(CPM + 1) of the up genes minus that of the down genes (CPM from the featureCounts total) | VERIFIED (`final_model_spec.txt`) |
| Genes | 9 genes, 8 up (KRT80, KRT23, GRIN2D, FUT1, ETV4, LARGE2, VSNL1, GTF2IRD1) and 1 down (OSR1); GENCODE v50 IDs | VERIFIED (`final_model_genes.tsv`) |
| Cutoff | **1.352325902**; call tumour if score > cutoff | VERIFIED (`final_model_spec.txt`) |
| Hash record | `lock/diagnostic_v1/sha256sums_final.txt` (sha256 `37657d53…6c`) covers 7 files: the 5 lock files plus `config/candidate_rules.yaml` and `config/reference_genes.yaml` | VERIFIED |
| Hash check | 7 of 7 OK on 2026-10-19 10:25 UTC | VERIFIED |
| Code freeze for validation | commit `495b4e6` (no tag); later commits up to `f4ab803` change docs only | VERIFIED |

The reference-gene list in `config/reference_genes.yaml` is covered by the hash record, but its formal approval is not recorded; it is
**not** described as locked here.

## 4. Tests and results

| Date | Test | Result | Status |
|--|--|--|--|
| 2026-10-15 | One-command diagnostic run on the 19 development samples (`scripts/run_diagnostic.sh`) | 48 of 48 regression checks PASS | VERIFIED (run log) |
| 2026-10-15 | `-entry DISCOVER` regression, 19 samples | 23 of 23 IDENTICAL or PASS; 1,390 candidates; same 9 genes and cutoff as the lock | VERIFIED (run log) |
| 2026-10-16 | TCGA-COAD technical validation of DISCOVER (499 public samples, 167 hold-out) | criteria D1 to D5 PASS; this tests the *method* on public data and fits a different 11-gene model, not the locked model | VERIFIED (run log) |
| 2026-10-16 | Format independence (renamed IDs, Tumor/Normal labels) | decisions identical to the lock; scores difference 0 | VERIFIED (`PHASE0_REPORT.md`) |
| 2026-10-19 | Pod baseline, recorded twice | decision checks 30 of 30 PASS each; run 2 vs run 1 byte-identical, 73 of 73 files | VERIFIED |
| 2026-10-19 | Regression after the folder-entry code | 73 of 73 files byte-identical to the baseline; 30 of 30 decision checks PASS | VERIFIED |
| 2026-10-19 | Input checks (bad files, labels, metadata) | 18 of 18 tests PASS | VERIFIED |
| 2026-10-19 | Folder entry, counts level, client-style IDs (RNI_001.., CohortA/B) | 16 of 16 checks against the lock PASS | VERIFIED |
| 2026-10-19 | Apply mode on the 19 development samples | locked scores identical (difference 0); development only | VERIFIED |
| 2026-10-19 | Resume of an interrupted run | 3 tasks reused from cache, 3 run, outputs identical | VERIFIED |
| 2026-10-19 | Smoke test, 1 real sample (SRR10967902), FASTQ to counts | **FAILED**: STAR killed for lack of memory at 29.60 GiB on a 29.8 GiB pod; no counts | VERIFIED |

Development performance of the locked model on its own 19 training samples: AUC 1.0, 19 of 19 calls correct (VERIFIED). This is
**development-only and optimistic** (genes and cutoff were chosen on these samples); it is not validation.

## 5. Known limits

- Research use only; not a clinical or diagnostic device.
- Human, paired-end, bulk RNA-seq only; single-end, multi-lane FASTQ and other species are not supported yet.
- Not built yet: Entrez or gene-symbol count input, a species guard, pathway analysis inside discover mode, a cohort outlier flag for
  very small apply cohorts.
- The locked model is valid only for colon tumour vs normal tissue, scored on featureCounts counts with GENCODE v50 IDs.
- STAR with the current index needs a pod with enough real memory: measured failure at 29.8 GiB (VERIFIED); recommended 64 GB RAM and 16
  CPUs, the configuration of the successful development run; the strict minimum was not measured (NOT VERIFIED).
- The development cohort is small (19 samples, 2 studies).

## 6. Current status

- **The independent 15-sample validation has NOT been run.**
- It is blocked by pod memory: the current pod has 8 CPUs and a 32,000,000,000-byte (29.8 GiB) memory limit (VERIFIED); STAR was killed
  during the smoke test (VERIFIED).
- Recommended pod: **64 GB RAM and 16 CPUs**, about 400 GB scratch space (NOT VERIFIED estimate for the scratch space).
- The 15 samples are ready: 7 normal and 8 tumour from two studies (SRP344867, SRP245232), 30 of 30 FASTQ files present, metadata
  checked (VERIFIED). None has been scored.

## 7. What happens next

1. Move to a pod with at least 64 GB RAM and 16 CPUs.
2. Pre-flight check: memory, CPUs, code version `495b4e6`, lock hashes 7 of 7.
3. Smoke test: one sample from FASTQ to counts and report, no scoring (about 40 minutes; NOT VERIFIED estimate).
4. With the owner's go: the 15-sample run in apply mode, locked model applied **once**, no retuning (about 5 to 8 hours; NOT VERIFIED
   estimate).
5. Report the pre-registered pass rule: at most 1 error per class, overall and within each study. AUC is reported, not the pass test.
6. Follow-up reports, report-only: rule 4 reproduction, and a reference-gene stability report (that script is not built yet).

A prognostic (survival) plan is drafted in `docs/PREREG_prognostic_DRAFT.md`; it is not approved and nothing has been run on real data.
