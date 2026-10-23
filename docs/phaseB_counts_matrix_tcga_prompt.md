# Task (Phase B): counts-matrix entry for the pipeline, then the TCGA-READ samples and the full 120-sample comparison with Hosseini 2023

Do NOT start until Phase A (docs/independent_validation_prompt.md) is finished and I say go.

Rules: same as Phase A. Do not modify main.nf, modules/, nextflow.config, bin/ or existing results without my approval of a written design. Write only to results/11_counts_matrix_tcga/ and new files under a new folder nextflow_counts_entry/. Compute every number from files. No 2>/dev/null. Deliver scripts by heredoc. Log in docs/run_log.md (Step 14) and docs/code_sources.md. Do not install anything on root /. If anything is unexpected, stop and report.

Step 0 (read-only, then STOP): read Table S1 of docs/reference/hosseini_2023_MOESM2_ESM.xlsx. Report exactly how the TCGA-READ samples (42 cancer / 43 normal) are identified (barcodes, UUIDs, or only counts) and what source they come from. Do not assume they are GDC primary tumour/normal samples; TCGA-READ normally has few normal samples, so check. If the sample identifiers or the source are not given, tell me and stop. Also check whether the pod can reach the GDC API (a small test only, no bulk download) and tell me what TCGA open-access data-use terms say about commercial use; report, do not decide.

Step 1 (design only, then STOP): propose a counts-matrix entry path: input = counts matrix + sample sheet (sample, group, study/source, baseline group named explicitly), an input check (baseline present, minimum samples per group, duplicate IDs, gene ID type; strip Ensembl version numbers and map to our GENCODE v50 IDs, report unmapped genes), then the same DESeq2 and pathway steps as the main pipeline with the same significance rule. It must be a separate entry that reuses the existing R scripts without editing them. Show me the design and the file list before writing code.

Step 2 (after approval): build it and test it first on our own 19-sample counts matrix (featureCounts output): results must equal the main run (same gene lists, same padj). Report any difference.

Step 3 (after approval): download the TCGA-READ counts identified in Step 0 and run the new entry on them. Report QC: sample counts per group, library size, PCA, outliers. TCGA counts come from a different upstream pipeline, so also run a combined analysis over all 120 samples with source/study as a covariate, and report per-source results separately.

Step 4: compare with the paper (Table S3, 5,037 genes) using the same columns as results/09_verification/step11_comparison_workbook (matched, opposite, missed, extra, recall, share of our hits also in the paper list) for the combined 120-sample analysis and for TCGA alone. Optional: if the metaSeq package installs without touching root /, run metaSeq across studies and report; if not, say so and skip.

Proposed pass rule (fixed before running; I may change it): for the combined 120-sample analysis, direction agreement >= 90% among genes in both lists. Report recall against the 5,037 list as information, no pass/fail.

Step 5: write results/11_counts_matrix_tcga/Phase_B_summary.md and a verified xlsx with a values-only copy, verify every number with a second independent script, and report.
