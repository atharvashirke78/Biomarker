# Session note 2026-10-12

## 1. Date and pod
- Date: 2026-10-12.
- This is a new pod, started after the pod stop on 9 Oct. The run was first interrupted by a pod stop on 8 Oct (about 17:16, recorded in run_log.md). The pod was stopped again around 10 Oct, after the last work at 12:16 on 9 Oct (exact time not in the log).
- Written read-only: nothing was run, installed or changed apart from this note and one line in docs/run_log.md. Versions below are copied from docs/tool_versions.tsv (recorded 2026-10-08) and were NOT re-checked on this pod.

## 2. Completed so far
- Diagnostic run: full 20-sample pipeline finished (all tasks done, 0 failed, 2026-10-09 09:28); docs/final_report_diagnostic.md written (19-sample DE: 10 normal, 9 tumour).
- T5_B (SRR13401190) conclusion: input is broken, not STAR. R2 is identical to R1 in the local files, in the ENA files (md5/bytes match ENA) and in a fresh SRA retrieval (prefetch + fasterq-dump, 8,885,228 reads per file = ENA read_count and Table S2). Sample FAIL in QC; no tumour replacement exists inside SRP301216.
- Verification steps 1-8 (results/09_verification/): per-sample QC, DESeq2 diagnostics, replication A vs B (BORDERLINE-weak), edgeR/limma-voom cross-check, Hosseini Table S3 direction/enrichment, sex check, VERIFICATION_SUMMARY.md, paired sensitivity + sample_map_with_patient.tsv.
- Per run_log.md, also logged on 2026-10-09: step8 extra questions (12:11) and verify_prompt2 remaining work incl. workbook move to docs/reference, step9_paper_comparison, leave-one-out T2_B (12:16).

## 3. In progress or missing when the pod stopped
- Part 2 checks (docs/verify_prompt2.md): verified, files exist (results/09_verification/step9_paper_comparison/, VERIFICATION_SUMMARY.md addendum 2).
- Hosseini workbook: verified, files exist; docs/reference/hosseini_2023_MOESM2_ESM.xlsx has md5 d260d2c65ef0293586f8b60792509de9 (matches) and code_sources.md has the line.
- step8 answers (results/09_verification/step8_extra_questions/, script scripts/step8_extra.R): verified, files exist.
- Step10 (2026-10-12): significant set redefined with the apeglm-shrunk |log2FC| >= 1 (results/09_verification/step10_shrunk_cutoff/, VERIFICATION_SUMMARY.md addendum 3). featureCounts 3,213 -> 2,798 genes, Salmon 3,281 -> 2,836 (Jaccard ~0.87, no gained genes, no direction flips); Table S3 direction/enrichment, the 18-gene table and ORA overlap with the paper are essentially unchanged.
- Not found in the docs: the printed answers 1-4 of the earlier state check (not available in this session).

## 4. Tools on this pod (base: /workspace/Biomarker_pipeline)
| tool | version (from tool_versions.tsv) | path |
|---|---|---|
| micromamba | 2.9.0 | env/bin/micromamba |
| Nextflow | (version not recorded) | env/bin/nextflow |
| Java | openjdk 17.0.18 | env/mamba/envs/nf/bin/java |
| R (r_dge) | 4.5.3; DESeq2 1.50.2, tximport 1.38.2, apeglm 1.32.0 | env/prefix/r_dge/bin/R |
| R (r_pathway) | 4.5.3; clusterProfiler 4.18.4, fgsea 1.36.2, msigdbr 26.1.1 | env/prefix/r_pathway/bin/R |
| R (edger env) | edgeR/limma (version not recorded) | /workspace/scratch_verify/edger/bin/R |
| samtools | 1.24 | env/prefix/samtools/bin/samtools |
| sra-tools | 3.4.1 | /workspace/scratch_T5B/sra/bin (prefetch, fasterq-dump) |
| fastp / FastQC / MultiQC | 1.4.0 / 0.13.0 / 1.35 | env/prefix/{fastp,fastqc,multiqc}/bin |
| STAR / featureCounts / Salmon | 2.7.11b / 2.1.1 / 1.10.3 | env/prefix/{star,subread,salmon}/bin |

```bash
B=/workspace/Biomarker_pipeline
export PATH=$B/env/bin:$B/env/mamba/envs/nf/bin:$PATH                 # micromamba, Nextflow, Java
export PATH=$B/env/prefix/r_dge/bin:$PATH                              # R (DESeq2 route); swap for r_pathway below
# export PATH=$B/env/prefix/r_pathway/bin:$PATH                        # R (pathway) - use instead of r_dge when needed
# export PATH=/workspace/scratch_verify/edger/bin:$PATH                # R (edgeR/limma)
export PATH=$B/env/prefix/samtools/bin:$PATH                           # samtools
export PATH=/workspace/scratch_T5B/sra/bin:$PATH                       # sra-tools
```
Note: both R installs provide `R`, so only one R directory should be first on PATH at a time.

## 5. Disk and memory
- /workspace (network FS): 2.1P total, 1.7P used, 412T available (81 % used). Biomarker_pipeline folder: 358 GB.
- Pod cgroup memory limit: 64,000,000,000 bytes (64 GB, ~59.6 GiB) (/sys/fs/cgroup/memory.max). cgroup usage at check time: ~1.6 GB. Host shows MemTotal ~755 GiB, which is not the pod limit.
- Nextflow config (from run_log): Salmon 28 GB x maxForks 2, STAR 33 GB x maxForks 1, inside a 58 GB executor budget.

## 6. Next steps (proposed order) and open decisions
1. Confirm part 2 / step8 / workbook status against the files (section 3) and re-check tool versions on this pod.
2. Hand VERIFICATION_SUMMARY.md and the step8 answers to the manager.
3. Re-run or finalise only what the manager decides after the points below.
Decisions waiting for the manager:
- KEGG licence (use of KEGG REST pathway sets / KEGG-derived results).
- Whether the 19-sample analysis (T5_B excluded) is accepted.

---
## Addendum (end of session, 2026-10-12): where we are before the pod stops

### Pod
- This pod: cgroup memory.max 16,000,000,000 (16 GB), nproc 4. Too small for the pipeline (needs memory.max >= 60 GB and 16 cores; Salmon 28 GB x 2, STAR 33 GB).
- `run_validation.sh` preflight (results/10_independent_validation/scripts/) was checked by hand and FAILS on this pod (memory and cores). Nothing was launched; no tmux run exists (only the Claude session).

### Validation run (15 samples, SRP344867 + SRP245232)
- Status: WAITING for a 64 GB / 16 vCPU pod. 15/15 downloads verified (integrity_check.tsv). Sample sheet validation_run/samplesheets/validation_15.csv ready.
- When the right pod is up: rerun the preflight; if it passes, launch in tmux with run_validation.sh, labelled EXPLORATORY (not the final validation) in the run log and results folder; after quantification, give the per-sample QC summary and stop (do not continue by itself). Note: the 19 + 15 combined validation on RC1 comes later (refactor Step 6).
- Not done: HISAT2. Feasibility was checked read-only (summary in the chat of 2026-10-12): the STAR route stays; no HISAT2 work now.

### Two-pipeline refactor (docs/two_pipeline_refactor_prompt.md)
- Step 0 (read only): DONE. Findings: memory guard in main.nf:86-90 only warns below 32 GB (silently drops to Salmon-only), no error and no --allow_salmon_only; no entry workflows; ~25 hard-coded thresholds in scripts/main.nf plus the ones already in nextflow.config; significance uses the UNSHRUNK log2FC; pathway code calls the KEGG web service (with MSigDB fallback).
- Git baseline: DONE. `git init` in /workspace/Biomarker_pipeline, .gitignore (results*, work*, data, test_data, logs, tmp, ref, .nextflow*, env/prefix|mamba|bin|nf_conda, large/binary files), commit 1e1f14b (83 files, code/config/docs only), tag `baseline-pre-refactor`. No existing file was changed. NOT in git: everything under results/ (incl. results/10_independent_validation/scripts/run_validation.sh) and docs/reference/*.xlsx. Files written after the commit (design doc, run_log, code_sources, these notes, .bak_* copies) are NOT committed yet.
- Step 1 (design document, no code): DONE, docs/two_pipeline_design.md. WAITING FOR OWNER APPROVAL. Steps 2-6 not started.

### Decisions made today (owner)
1. Aligner stays STAR (+ featureCounts) with Salmon/tximport; no HISAT2 now.
2. One pipeline, two stages (QUANT, ANALYSIS); no -entry runs both.
3. Significant genes: padj < 0.05 and |apeglm-shrunk log2FC| >= 1; unshrunk value also reported.
4. Pathway background = genes tested by DESeq2 after the low-count filter.
5. Sample sheet gets count_source, baseline, batch; group names and baseline come from the sheet.
6. cohort_qc pre-filter and DESeq2 pre-filter stay as two named parameters with current values (10 counts in min(3,n) samples; 10 counts in 5 samples).
7. Pathway gene sets from local versioned files recorded in the manifest (name, version, md5); no live web calls.
8. Validation-only scripts (check_vs_paper.R / .nf) move out of the product path.
9. KEGG stays on but switchable (--run_kegg); licence note in report and docs.
10. Below 32 GB (or unreadable limit) main.nf must stop with an error unless --allow_salmon_only true; method used recorded in manifest and report header.
11. Regression test rules fixed: ANALYSIS-only on the existing 19-sample matrix; QUANT-only on the 4-sample mini test (needs 64 GB pod and approval); no 19-sample STAR/Salmon re-run until the combined 19 + 15 run on RC1.

### Open decisions (marked TBD in docs/two_pipeline_design.md section 10)
min samples per group; GSEA rank metric (Wald stat proposed); GO simplify cutoff; Reactome/KEGG release and KEGG GMT source; max unmapped gene % check; keep legacy 4-column sample sheet or not; >2 groups out of scope for v1; concordance.R placement; whole-index md5 vs size; whether run_validation.sh should be copied into scripts/ and versioned. Also still open from before: KEGG commercial licence, and acceptance of the 19-sample analysis (T5_B excluded).

### Next steps
1. Owner approves (or changes) docs/two_pipeline_design.md and answers the TBDs; then Step 2 (implement entry workflows, ANALYSIS validation, manifest writer) and Step 3 (central thresholds file, STOP to show the list).
2. Start a 64 GB / 16 vCPU pod; run the preflight; then the EXPLORATORY 15-sample validation run as above. Do not edit pipeline files while a run is in progress.
