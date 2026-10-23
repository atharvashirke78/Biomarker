# Two-pipeline refactor, config centralisation, pathway setup (STAR route stays)

## Decisions already made (do not reopen)
- Aligner stays STAR (+ Salmon, featureCounts). No HISAT2 work now.
- One pipeline, two internal stages in the same repo:
  - QUANT: raw FASTQ -> fastp QC -> STAR/featureCounts + Salmon/tximport -> per-sample QC -> counts matrix + QC table + MultiQC + manifest.
  - ANALYSIS: counts matrix + sample sheet -> input validation -> normalisation -> DESeq2 DEG -> pathways -> diagnostic report. Prognostic/predictive branches come later; leave a clean place for a selector (--analysis).
- Running without -entry must still run QUANT then ANALYSIS, same as today.
- Significant genes: padj < 0.05 and |shrunk log2FC| >= 1 (apeglm). Also report the unshrunk value.
- Pathways: GO (BP, MF, CC; simplified), KEGG, Reactome, optional MSigDB Hallmark. For each: over-representation (background = all expressed genes) AND GSEA on the ranked list (shrunk log2FC or Wald stat). BH correction per database, databases reported separately, never pooled.
- KEGG stays ON for now but must be a switchable module (--run_kegg true/false). Licence for commercial use is unresolved: print a clear note in the report and in docs.
- If not already done: below 32 GB memory (or unreadable limit) main.nf must stop with an error, unless --allow_salmon_only true is passed. Record the method actually used (STAR+featureCounts or Salmon-only) in run metadata and report header.

## Rules
- Do not launch Salmon or STAR jobs without my approval. Code editing and -stub-run are fine on any pod; the mini test and full run need the 64 GB / 16 vCPU pod (run_validation.sh preflight shows what the pod has).
- Do not delete or overwrite anything under results/. Keep .bak copies of every file you edit.
- No 2>/dev/null on first runs. fastp QC must be in every results folder. Per-sample folders named by SRR.
- Log every change in docs/run_log.md and docs/code_sources.md (code provenance).
- Do not edit pipeline files while a pipeline run is in progress.
- Stop at each STOP and wait for me.

## Step 0 (read only)
Read the session notes md and docs/ in /workspace/Biomarker_pipeline, nextflow/main.nf, nextflow.config, modules/ and bin/. Report: current state, whether the STAR memory hard-error is already in, and every hard-coded threshold you find (file and line). STOP.

## Step 1 design document, no code
Write docs/two_pipeline_design.md with the file contract between QUANT and ANALYSIS:
- counts_matrix.tsv (genes x samples, raw integer featureCounts), tximport gene counts and TPM, sample_sheet.tsv (sample_id, group, study, batch, count_source, baseline flag), qc_table.tsv (pass/fail per sample with reason), manifest.json (pipeline version, tool versions, all resolved parameters, md5 of inputs and references).
- ANALYSIS input checks that stop the run with a clear message: baseline group missing, fewer than min samples per group (parameter, value to be agreed, mark TBD), non-integer counts, sample sheet and matrix mismatch, count_source not recorded.
- Entry workflows: -entry QUANT, -entry ANALYSIS, default both.
STOP for approval.

## Step 2 implement
Named entry workflows; ANALYSIS validation; manifest writer.

## Step 3 central config
Create one file (e.g. conf/thresholds.config) holding every threshold (padj, log2FC, low-count filter, QC warn/fail limits, min samples per group, pathway cutoffs). Nothing hard-coded in scripts. Write the resolved values into manifest.json. STOP and show me the list.

## Step 4 DE and pathways
Shrunk-LFC cutoff and the pathway design above. Add a short Methods text for each database with tool versions.

## Step 5 regression test (pass rules, fixed now)
Do NOT re-run the 19 samples through STAR/Salmon for the regression test.
- (a) Run ANALYSIS only, starting from the existing 19-sample counts matrix, and compare with the existing DESeq2 results. Pass: identical except the intended changes (shrunk cutoff, new gene sets), which must be listed gene by gene.
- (b) Run QUANT only on the 4-sample mini test (needs the 64 GB pod, after my approval) and compare its counts with the old pipeline's counts for the same samples. Pass: exact integer match.
- (c) The full 19-sample re-run happens once, inside the combined 19 + 15 validation run on RC1.
Report pass/fail per rule. STOP.

## Step 6 release candidate
If Step 5 passes: tag RC1 (git tag, or a VERSION file if no git), freeze env lock files, record md5 of genome, GTF v50 and indexes, write docs/release_RC1.md. Do not start the combined 19 + 15-sample validation run; wait for my go.
