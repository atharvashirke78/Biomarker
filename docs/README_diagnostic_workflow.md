# Diagnostic workflow: how to run it (rc1-prep, lock-diagnostic-v1)

Colorectal cancer, tumour vs normal. The model is locked at tag `lock-diagnostic-v1` (9 genes, balanced up/down log2-CPM score, cutoff 1.352325902).

- Nothing in this workflow refits or changes the model.
- Development results (19 samples from 11 patients) are optimistic.
- No independent validation has been run yet.

Run every command from `/workspace/Biomarker_pipeline`. Output folders must not exist: the scripts never overwrite.

## 1. Development workflow (runs on the current pod)

Starts from the existing quantification of the 19 development samples (`results/04_counts/counts_featurecounts.tsv`, `results_rc1prep_19/sample_sheet_rc1prep.tsv`). It doesn't need FASTQ, STAR or Salmon. Takes about 10 min on 8 cores and less than 2 GB RAM.

```bash
scripts/run_diagnostic.sh results_rc1prep_19/diagnostic_run_<YYYYMMDD_HHMM>
```

What it runs:
1. `nextflow -entry ANALYSIS`: input checks, then DESeq2 with `~ batch + group`, then GO BP and Hallmark ORA and GSEA. KEGG runs only with `--run_kegg true` (default `false` since commit `8db47d6`).
2. `scripts/apply_locked_model.R`: applies the locked score to the 19 samples, labelled "development, optimistic, not validation".
3. `scripts/rule4_reproduction.py`: a self-consistency test, which must reproduce all 1,390 candidates.
4. `tests/regression/check_19_development.py`: compares the outputs with the documented expected values, then confirms the new DESeq2 table is byte-identical to `results_rc1prep_19/deseq2_results_all.tsv`.

Single steps:
```bash
# locked model on any featureCounts matrix (versioned GENCODE v50 gene IDs)
env/prefix/r_classifier/bin/Rscript scripts/apply_locked_model.R --counts <counts.tsv> --cohort development|validation \
    [--sample-sheet <contract sample_sheet.tsv>] [--qc-table <qc_table.tsv>] [--evaluate] --out-dir <new dir>
# regression checks on the existing 19-sample outputs
env/prefix/pyqc/bin/python3 -I tests/regression/check_19_development.py --out <summary.tsv>
```

Safeguards in `apply_locked_model.R`:
- It checks the lock-file sha256 against `lock/diagnostic_v1/sha256sums_final.txt` before scoring.
- It reads the genes, signs and cutoff from the locked files.
- It stops if any model gene is missing from the matrix.
- It drops QC FAIL samples automatically and lists them.
- `--evaluate` on validation data requires the QC table.

## 2. Independent validation (15 held-out samples): not yet runnable here

The FASTQs are present (`validation_run/data/`, 30 files, 96.1 GB, sizes match ENA). The pod needs at least 60 GB cgroup memory and 16 cores. Check the current pod with `scripts/run_validation.sh preflight`; it failed on 2026-10-12 (16 GB, 4 cores) and on 2026-10-15 (32 GB, 8 cores).

```bash
scripts/run_validation.sh preflight      # checks only; PASS is required before anything starts
scripts/run_validation.sh quant          # -entry QUANT in the background; prints the tail -f and resume commands
tail -f results_rc1prep_15quant/quant_run.log
scripts/run_validation.sh resume         # if the run stops
scripts/run_validation.sh downstream     # after QUANT: ANALYSIS, locked model (applied once), rule 4
```

The memory guard is never bypassed. Validation AUC, sensitivity and specificity come only from `downstream`, which runs after all 15 samples are quantified with verified labels. The locked pass rule allows at most 1 error per class, overall and within each study.

## 3. Discovery steps (locked; re-run only to reproduce, never to change the model)

| Step | Command | Note |
|--|--|--|
| Candidate rules 1-3b | `env/prefix/r_dge/bin/Rscript scripts/candidate_selection.R results/04_counts results_rc1prep_19 <new dir>` | Guarded to the 19-sample cohort (rule 3b = 18 of 19); about 3 min |
| Reactome ORA/GSEA | `R_USER_CACHE_DIR=$PWD/env/msigdbr_cache env/prefix/r_pathway/bin/Rscript scripts/downstream/reactome_analysis.R results_rc1prep_19 <new dir>` | Thresholds hard-coded to equal `conf/thresholds.config` |
| LASSO / lock fit | `scripts/downstream/lasso_*.R`, `results_rc1prep_19/lock_proposal_v2/fit_final_model_v2.R` | Tied to the 19 samples; the model is locked, so these are not part of the workflow |
