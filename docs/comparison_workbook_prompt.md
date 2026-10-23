# Task: build the paper-vs-pipeline comparison workbook (READ-ONLY on all existing results)

Goal: one Excel workbook that Sakshi can show to her manager (Jaspal Patil). For every comparison against Hosseini et al. 2023 (docs/reference/hosseini_2023_MOESM2_ESM.xlsx) it must show: what the paper reports, what our pipeline reports, how many we matched, how many we missed (and WHY), and how many are extra in ours.

Rules:
- Do not modify any existing result or pipeline file. Write only into results/09_verification/step11_comparison_workbook/.
- Compute every number from the real files. Do not type numbers from memory. If a number differs from the checks listed in the last section, STOP and report it, do not silently fix it.
- Do not suppress errors with 2>/dev/null. Do not install anything on root /, which has about 4.8 GB free. If openpyxl is missing: pip install --target /workspace/scratch_verify/pylib openpyxl and add that folder to PYTHONPATH.
- Write the build script to results/09_verification/step11_comparison_workbook/build_comparison.py so it can be rerun.
- Append a line to docs/run_log.md (Step 12, comparison workbook) and a line to docs/code_sources.md for the new script (own code; libraries: pandas, openpyxl).

Inputs (find exact paths with ls first and show me what you used):
- Paper: Table S2 (read counts per run), S3 (DEGs up/down with probability, no fold change), S4/S5 (GO), S6 (KEGG), the 18-gene table used earlier in step9.
- Ours: DESeq2 results from featureCounts and from Salmon (final 19-sample run; significant = padj < 0.05 and |log2FC| >= 1, unshrunk, as in deseq2.R), our GO and KEGG ORA results, step9_paper_comparison files.
Use featureCounts-based DESeq2 as the primary column and Salmon as a secondary column.

Workbook name: Hosseini2023_vs_Pipeline_comparison.xlsx. Sheets:

1. Summary. One table with the columns: Item | Paper | Our pipeline | Matched | Missed (paper only) | Extra (ours only) | Recall % | Our hits also in paper list % | Note. Rows: samples used (paper 120 total; ours 19 = 9 tumour + 10 normal, 20 taken, T5_B excluded because R2 is a copy of R1); read-count agreement with S2 (runs matched out of 20); down genes; up genes; all DEGs; GO terms; KEGG down pathways; KEGG up pathways; the 18-gene table (significant / expected direction). Use Excel formulas for recall and for 'Our hits also in paper list %' (matched / ours) so they are live. The column is NOT called precision: the paper used metaSeq on 120 samples, so ours-only genes are not necessarily false positives; say this in the Note column of the gene rows and in the 'How to read this' box. Below the table add a short plain-English "How to read this" box: agreement with the paper is not independent validation, because our 20 samples are part of the paper's data; low up-gene recall is expected from the 19 vs 120 sample size, and that explanation is stated as likely, not proven.

2. Gene_overlap. One row per gene in the union of (paper DEGs, our significant genes): gene | paper direction | our direction | our log2FC (unshrunk) | our padj | our baseMean | category. Categories: Matched (same direction); Opposite direction; Paper only: tested, padj >= 0.05; Paper only: padj < 0.05 but |log2FC| < 1; Paper only: padj NA (Cook's / low count); Paper only: not in our DESeq2 table (filtered or not annotated); Ours only (not in paper list). Add counts of each category at the top, with autofilter and frozen header.

3. Missed_breakdown. A small table of the category counts for up and for down, with percentages of the paper list. This answers "how much did we miss and why". Add one bar chart (openpyxl) of matched vs each missed category.

4. Pathways_GO_KEGG. Paper term | direction | in our results? | our adjusted p | note. Matching on KEGG pathway name for KEGG; for GO match by GO ID. List clearly the terms found only in the paper (include spliceosome, RNA transport, mRNA surveillance if confirmed) and terms only in ours. Note that the KEGG release difference is a likely explanation and is not verified.

5. Readcounts_S2. Run | paper total reads | our total reads | match (TRUE/FALSE); 20 rows; T5_B marked as excluded from DE analysis.

6. Gene18. The 18-gene table: gene | paper direction | our direction | our padj | significant? | caveat (CXCL8, KRT6A, KRT6B, KRT17 direction-only because padj NA; COL3A1 borderline; APOA4 untestable).

7. Method_notes. Pipeline versions (fastp 1.4.0, STAR 2.7.11b, Salmon 1.10.3, featureCounts, DESeq2 design ~ study + group), reference (hg38 chr1-22,X,Y,M; GENCODE v50), the significance rule, the sample list, and the file paths used for every sheet.

Formatting: header row bold with fill, freeze panes, column widths set, numbers formatted (padj in scientific notation, percentages as 0.0%), no merged cells in data tables, landscape page setup, sheet names exactly as above.

Verification (all required before reporting):
a) Reopen the saved workbook with openpyxl and with pandas; confirm there are no formula error values and all 7 sheets exist.
b) Recompute the Summary numbers independently from the raw files in a second small script and confirm they equal the workbook values.
c) These known values must be reproduced (corrected 2026-10-12; computed from the files, featureCounts DESeq2, significant = padj < 0.05 and |log2FC| >= 1 unshrunk): S2 read counts matched 20/20; paper lists 285 down and 4,752 up; our significant genes 2,060 down and 1,153 up (gene-id level; 1,152 up at gene-symbol level because two Ensembl IDs share one symbol); of these 223 down and 865 up are in the paper's lists (223 in the down list, same direction; 865 in the up list, of which 855 with the same direction and 10 opposite); enrichment 7.2x down and 3.9x up; KEGG down 11/23 and up 5/8 shared. Stop only if the files disagree with these values, and tell me which and by how much. (The earlier wording 'our 223 significant down' / '855 of our 865' mixed up the counts and is replaced by the above.)
d) Check that the sum of the categories in Gene_overlap equals the union size, and that Matched + Opposite + Paper-only categories for each direction add up to the paper list size.
e) Show me the Summary sheet as a plain-text table, and the category counts.

Finally print the file path and size, and the one-line command to send it to the laptop: runpodctl send <path>.
