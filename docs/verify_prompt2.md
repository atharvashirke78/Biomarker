# Verification part 2: compare with Hosseini et al. 2023 supplement and finish the checks

Rules: read-only on all pipeline code and earlier results. Write only into results/09_verification/ (scripts in results/09_verification/scripts/). The only other allowed changes are step 0. Do NOT rerun the pipeline and do NOT change any cutoff or setting to improve agreement. Cutoffs stay padj < 0.05. Do not suppress errors with 2>/dev/null.

## Step 0: place the reference file
The Hosseini supplementary workbook was uploaded to /workspace/Biomarker_pipeline/41598_2023_44459_MOESM2_ESM.xlsx. Check md5sum is d260d2c65ef0293586f8b60792509de9 and size 1291527 bytes. If it matches: mkdir -p docs/reference, then mv -n it to docs/reference/hosseini_2023_MOESM2_ESM.xlsx, and add one line to docs/code_sources.md (source, file name, md5, date). If it does not match, stop and report. From now on use this file instead of the copy downloaded from Springer.

## What the workbook contains (verified)
- Table S1: the paper integrated 4 SRA studies (SRP219837 5v5, SRP301216 5v5, SRP344867 5v4, SRP245232 3v3) plus TCGA READ (42 cancer, 43 normal); 120 samples in total. We reproduce only 2 SRA studies (20 samples), so only direction-level agreement is expected. State this in the summary.
- Table S2: per-run Total Reads for each SRR. SRP301216 lists SRR13401190 (8,885,228 reads) as a cancer sample, so the paper used T5_B.
- Table S3: Upregulated genes (cols A-C, 4,752) and Downregulated genes (cols E-G, 285). It has a metaSeq probability and p-value, NO fold change. So compare direction only (listed up / listed down / not listed).
- Table S4 and S5: GO terms for down and up genes (BP, MF, CC). Table S6: KEGG for down and up genes.

## Step 1: read-count integrity against Table S2
For all runs of SRP301216 and SRP219837 in sheet Table S2, compare our fastp raw read counts (read pairs) per SRR with S2 Total Reads. Print a table: SRR, our reads, S2 reads, match yes/no. Include SRR13401190 (S2 says 8,885,228). Say whether any run differs and by how much.

## Step 2: direction comparison with Table S3
Use our DESeq2 featureCounts result (results/06_dge/deseq2_results_all.tsv) and also the Salmon result. Map Ensembl IDs to gene symbols with the GENCODE v50 annotation already used by the pipeline; report how many S3 symbols could not be matched (symbol aliases such as H1-3, KCNQ1OT1 etc.) and do not hand-fix them silently.
Report, as text tables:
(a) Universe = genes with non-NA padj in our DESeq2. Our significant up genes (padj < 0.05, log2FC > 0): how many are in S3 up, in S3 down, not listed. Same for our significant down genes. Direction agreement among those that are listed.
(b) Hypergeometric enrichment: our significant up vs S3 up, our significant down vs S3 down, with the background restricted to genes present in both our tested universe and S3 symbols. Give overlap, expected, odds ratio, p.
(c) Concordance over all S3-listed genes present in our results: percent with sign(log2FC) matching the list, and the same restricted to our padj < 0.05 and to the top 100 S3 genes by metaSeq in each list.
(d) The 18-gene table as plain text: gene, our log2FC, our padj, our direction, S3 list (up/down/not listed), same direction yes/no.

## Step 3: pathway comparison with S4, S5, S6
Compare our ORA results (results/07_pathway/ora, KEGG and GO BP) with S6 (KEGG up and down) and S4/S5 (GO BP). List terms found in both (match by name), terms only in the paper, and terms only in ours, with the top 10 of each. Note that the KEGG releases differ.

## Step 4: sensitivity analyses (new files only)
(a) Paired analysis for study A: patient IDs come from the ENA sample titles (e.g. COAD-T004). The labels T1_A, N1_A are NOT patient numbers. Use design ~ patient + group, ref normal, on the real pairs. Check ENA titles for study B too and say whether pairs exist there. Compare with the main unpaired result: number significant, Jaccard of significant sets, Spearman of log2FC, the 18-gene table.
(b) Leave-one-out for SRR13401202 (T2_B, flagged by cohort QC with robust z -11 in a cell of 4): rerun DESeq2 without it (same design ~ study + group, ref normal, apeglm) and compare with the main result with the same four metrics plus direction agreement with Table S3.
(c) Write results/09_verification/sample_map_with_patient.tsv as a NEW file with a patient_id column. Do not edit sample_map.tsv.

## Step 5: T5_B
Report the result of the SRA retrieval for SRR13401190: is R2 identical to R1, and what is the read count compared with ENA and with S2 (8,885,228)? Do not rerun anything.

## Step 6: report
Update results/09_verification/VERIFICATION_SUMMARY.md with all of the above, with a short "Expected vs observed" note (paper = 4 SRA studies + TCGA; we have 2 SRA studies). Then print the new sections as plain text in your reply, in chunks of about 50 lines, so the output does not collapse. Also print the replication numbers (Spearman, Pearson, n genes, % same direction, study A vs study B), the edgeR/limma cross-check numbers, and list the PNGs copied to figures_for_download/.
