# Task: independent validation of the diagnostic branch on two unused studies (SRP344867, SRP245232)

Purpose: test the pipeline on data it has never seen. The paper (Hosseini 2023) lists SRP344867 (5 cancer, 4 normal) and SRP245232 (3 cancer, 3 normal). Our 19-sample analysis used two other studies.

Rules:
- Do NOT modify main.nf, modules/, nextflow.config, bin/ or any existing result. If the pipeline cannot take a different sample sheet and output folder through parameters, STOP and tell me; do not edit it.
- Run from a separate launch folder (/workspace/Biomarker_pipeline/validation_run/) with its own work folder (/workspace/validation_work). Do NOT use -resume of the main session and do not touch .nextflow/history of the main run.
- All results go to results/10_independent_validation/<study>/ with per-sample folders named by SRR and the fastp QC present in each.
- Do not suppress errors with 2>/dev/null. Deliver scripts by heredoc. Compute every number from files, never from memory. If anything is unexpected, stop and report; do not fix silently.
- Log in docs/run_log.md (Step 13) and docs/code_sources.md for every new script. Do not install anything on root /.

Step 0 (read-only): check how the main pipeline takes its sample sheet (columns, format) and outdir, and tell me whether it can run on a different sample sheet without edits.

Step 1 (metadata, no download yet): get run accessions and metadata for both projects from ENA (sample titles, tissue, condition, library layout, read length, read counts, sizes). Show a table per study. Check: cancer/normal counts equal the paper's (5/4 and 3/3); whether samples are paired by patient; tissue type and library prep; whether runs are listed in paper Table S2 and, if so, whether read counts match. If a study is not primary tissue RNA-seq or the labels are unclear, say so. STOP here and wait for my approval before downloading.

Step 2 (after approval): download, run the identical R1/R2 check and fastp QC, then run the pipeline on each study as a separate cohort (design ~ group, reference normal). Apply the same QC rules (FAIL excludes a sample). Report the per-sample QC status.

Step 3 (comparison, per study): (a) direction agreement with paper Table S3 among genes significant in the new study and in the paper's list; (b) overlap and enrichment against the paper's lists; (c) Spearman correlation of log2FC with our 19-sample DESeq2 results over all shared genes, and direction agreement among genes significant in both; (d) the 18-gene table; (e) KEGG/GO overlap with paper S4-S6. Use the same significance rule as the main run (padj < 0.05 and |log2FC| >= 1).

Proposed pass rules (written before the results; I may change them): for SRP344867, direction agreement >= 90% among genes significant in both the new study and the paper list. For SRP245232 (3 vs 3, low power) report only, no pass/fail. State at the top of the report that these were fixed before the run.

Step 4: write results/10_independent_validation/Independent_validation_summary.md and an xlsx (Summary, Gene_overlap per study, Pathways, Samples_QC, Method_notes) with a values-only copy. Verify every number with a second independent script and report counts.
