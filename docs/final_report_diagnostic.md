# Diagnostic full run (20 samples) - final report

Pipeline: `nextflow/main.nf` (profile `conda`), session `d6ae67d4-b110-462a-84f6-78b1e3219408` (`docs/full_run_session_id.txt`).
Run log with all events and decisions: `docs/run_log.md`. Tool versions: `docs/tool_versions.tsv`.

## 1. Samples used: 19 of 20

20 samples (10 tumour, 10 normal) from two studies in `samplesheets/diagnostic_20.csv`: SRP219837 (study A, 5 T + 5 N) and SRP301216 (study B, 5 T + 5 N).

**Excluded: T5_B (SRR13401190, SRP301216, tumour) - per-sample QC status FAIL** (`results/04_counts/sample_qc.tsv`).

| value | T5_B | rule |
|---|---|---|
| STAR uniquely mapped | **0.02 %** | FAIL below 50 (WARN below 70) - this triggered the FAIL |
| featureCounts assigned | 0.06 % | WARN below 60 |
| Salmon mapped | 60.68 % | WARN below 70, FAIL below 50 (not failed) |
| fastp Q30 after filtering | 89.23 % | - |
| fastp reads passing filter | 98.17 % | - |

Why: the **input data is broken, not STAR**. In the raw files `data/tumor_sample/SRR13401190_R1.fastq.gz` and `_R2.fastq.gz` the read-2 sequences are identical to the read-1 sequences (all 20,000 reads compared at the start; the md5 of all sequence lines of the full files is the same for R1 and R2). The other 19 samples have 0 of 5,000 identical pairs. Consequences seen in the tools:
- STAR: 99.02 % of read pairs "unmapped: too short" (a pair of two identical mates cannot align as a proper pair), 2,140 unique pairs.
- Salmon: 5,293,319 fragments mapped, but only 175 concordant; 5.58 M orphan/inconsistent (library counts SF 2.83 M / SR 2.75 M, essentially no IU).
- fastp: identical adapter/duplicate statistics for read 1 and read 2 (duplication 63.6 %).
The pipeline's own VERIFY_FASTQ (read counts R1 = R2 = 8,885,228, gzip integrity) cannot detect this. Not fixed here: the run accession should be re-downloaded and re-checked before T5_B is used (study B tumour then has 4 samples).

Per-sample status of all 20: PASS 14, WARN 5 (T3_A, N1_A, N4_A, N5_A, N5_B), FAIL 1 (T5_B); no WARN sample was removed. All cutoffs unchanged.

## 2. Settings (nextflow/nextflow.config)

- Reference: GRCh38 primary chromosomes, GENCODE v50 (GTF + transcripts), chrY PAR genes removed from the gene set; STAR index sjdbOverhang 75 (from T1_A; study B reads are 101 bp).
- fastp: Q20, 40 % unqualified limit, cut_tail window 4 / mean Q20, min length 36, max 5 N bases. STAR: `outFilterMultimapNmax` 20. featureCounts: `-Q 10`, strandedness from the Salmon library type.
- Salmon: decoy-aware index, `--gcBias --seqBias`, validateMappings, libType A.
- DESeq2: design `~ study + group`, reference `normal`, pre-filter >= 10 counts in >= 5 samples, significance padj < 0.05 and |log2FC| >= 1 (unshrunk). Two routes: featureCounts matrix and Salmon tximport.
- Pathways: ORA (clusterProfiler, KEGG REST) and fgsea (msigdbr 2026.1.Hs: Hallmark, C2 CP:KEGG_LEGACY, GO BP), set size 10-500, q 0.05, 10,000 permutations.
- Resources: executor local 16 CPUs / 58 GB (cgroup 64 GB). STAR 33 GB, `maxForks 1`. Salmon changed on 2026-10-09 to 28 GB fixed, 8 CPUs, `maxForks 2` (2 x 28 = 56 GB); measured Salmon peak_rss 26.2-26.7 GB, STAR 29.0-29.1 GB.

## 3. Tool versions (docs/tool_versions.tsv)

Nextflow 24.10.5, Java 17.0.18, micromamba 2.9.0, fastp 1.4.0, FastQC 0.13.0 (not run, `run_fastqc = false`), MultiQC 1.35, samtools 1.24, STAR 2.7.11b, featureCounts 2.1.1, Salmon 1.10.3, Python 3.12.15 + pandas 2.3.3, R 4.5.3 with DESeq2 1.50.2, tximport 1.38.2, apeglm 1.32.0, clusterProfiler 4.18.4, fgsea 1.36.2, msigdbr 26.1.1, org.Hs.eg.db 3.22.0, enrichplot 1.30.5.

## 4. Run times and the pod stop

- Final (resumed) session: started 2026-10-09 06:32, completed 09:08:49, **duration 2h 36m 34s**, **118 tasks succeeded, 62 cached**, CPU hours 46.0 (25.5 % cached), no failed task (`logs/nf_full_20261009_063206.log`, trace `logs/nf_trace_20261009_063206.tsv`).
- The first part of the run (2026-10-08, 15:1x-17:16) is in `logs/nf_full_part1.log`, `logs/old_20261009_0627/` and `logs/nf_*_part1.*`.
- **Pod stop:** the container was restarted at about 17:16-17:17 on 2026-10-08 (external; no OOM/exit 137 in the trace or log; the tmux session `nf_diag` was lost). Lost: the two running tasks (STAR_ALIGN T3_A, FEATURECOUNTS T3_B) and the contents of `.nextflow/history` (0 bytes). `work/` and `results/` were intact.
- **Resume:** `-resume d6ae67d4-...` failed first ("Can't find a run with the specified id") because `history` was empty; after backing up `.nextflow` and the logs, one tab-separated line was written to `.nextflow/history` and the resume worked: FASTP 20, STAR 3, Salmon 6, featureCounts 2 were cached.
- Config change applied at the resume: Salmon 28 GB, 8 CPUs, `maxForks 2` (13 Salmon pairs overlapped in time, no OOM, no failure).

## 5. Differential expression (19 samples: 10 normal, 9 tumour; `results/06_dge`)

| route | tested | padj < 0.05 | up | down | up + down (padj < 0.05 and |log2FC| >= 1) |
|---|---|---|---|---|---|
| featureCounts | 17,030 | 4,398 | 1,153 | 2,060 | 3,213 |
| Salmon (tximport) | 18,112 | 4,009 | 1,123 | 2,158 | 3,281 |

Concordance featureCounts vs Salmon (`concordance_featurecounts_vs_salmon.tsv`): 16,870 genes in both; significant overlap 2,554 of union 3,940, **Jaccard 0.65**; log2FC Pearson **r = 0.95** over all common genes (Spearman 0.93), r = 0.97 for genes significant in either route.

## 6. Check against the paper (after the fact; no cutoff changed)

`results/06_dge/check_vs_paper_*.tsv`:
- Genes: **direction matches 17 of 18**; direction and padj < 0.05: **12 of 18**; 1 not tested.
- KEGG pathways: **10 of 12** match (ORA or GSEA significant in the expected direction). The two that do not are RNA transport and mRNA surveillance (expected up; not significant in ORA, and not available in GSEA, see below).
- GSEA KEGG: 7 of the 12 pathways show NA for GSEA (mineral absorption, fat digestion and absorption, pancreatic secretion, bile secretion, RNA transport, mRNA surveillance, IL-17 signaling). **Cause: the gene set is missing, not a name-matching error.** fgsea uses msigdbr `C2 / CP:KEGG_LEGACY` (186 sets; 184 passed the 10-500 size filter), and that collection has none of these 7 pathways (the only partial name hits are KEGG_PANCREATIC_CANCER and KEGG_PRIMARY_BILE_ACID_BIOSYNTHESIS, which are different pathways); the other KEGG collection in msigdbr, CP:KEGG_MEDICUS (658 sets), has no matching names either. ORA uses the current KEGG REST sets and does contain 5 of them (mineral absorption, fat digestion and absorption, pancreatic secretion and bile secretion down; IL-17 signaling up), so those 5 are reported as matches through ORA. RNA transport and mRNA surveillance are not enriched in the ORA up list. The check script (`nextflow/bin/check_vs_paper.R`, `norm()` + fixed-string match) finds the pathways that do exist (e.g. PPAR signaling, cell cycle, spliceosome, DNA replication).

## 7. Cohort QC / PCA (`results/05_cohort_qc`)

- PC1 (48.0 % of variance) = **study**, R2 0.84 (group R2 0.05). PC2 (14.3 %) = **group**, R2 0.55 (study R2 0.10). PC3 10.8 %, PC4 7.4 %.
- Outlier flag (flag only, nothing removed) changed on 2026-10-09: computed within each study x group cell (robust z > 3 MAD on PC1 or PC2) instead of within group across both studies. Before: 6 samples of study B flagged by the batch effect and the JSON showed none (it read a column that did not exist). Now: TSV and JSON agree, **1 flagged: SRR13401202 (T2_B)**, PC1 robust z = -11 (cell of only 4 samples, so the z is unstable; treat as "look at it", not as an exclusion). Previous files: `results/05_cohort_qc/old_all_samples_flag_20261009_0918/`. Script backup: `nextflow/bin/cohort_qc.R.bak_pre_within_study`. The change was applied by re-running only `cohort_qc.R` on the original inputs (PC1/PC2/library sizes identical), not by re-running the pipeline; a later `-resume` will not recompute it automatically.

## 8. Results folder map (`results/`)

| folder | content |
|---|---|
| `01_fastp/<SRR>/` | fastp.html / fastp.json per sample; `qc/` copies of all fastp reports + `multiqc_report.html` |
| `02_alignment/<SRR>/` | `<SRR>.sorted.bam(.bai)`, `.flagstat`, `.stats`, STAR `Log.final.out`; `qc/` |
| `03_quant_salmon/<SRR>/` | Salmon quant directory (`quant.sf`, `aux_info/`, `lib_format_counts.json`); `qc/` |
| `04_counts/<SRR>/` | featureCounts table and summary per sample; cohort `sample_qc.tsv`, count matrices; `qc/` |
| `05_cohort_qc/` | PCA (group/study), correlation heatmap, library sizes, `cohort_qc_flags.tsv`, `cohort_qc_summary.json`; `qc/` |
| `06_dge/` | DESeq2 results/normalised counts/significant up and down/volcano/MA/heatmaps for both routes, `concordance_featurecounts_vs_salmon.tsv`, `check_vs_paper_*`; `qc/` |
| `07_pathway/ora/`, `07_pathway/gsea/` | ORA (KEGG REST, GO BP) and fgsea (Hallmark, KEGG legacy, GO BP) tables and plots; `qc/` |
| `08_multiqc/` | cohort MultiQC report |
| `sample_map.tsv`, `RESULTS_INDEX.md` | SRR / label / study / group map; index of the folders |

Size: results 35 GB, work 138 GB (not deleted).
