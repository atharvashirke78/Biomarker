# Two-stage pipeline design (QUANT / ANALYSIS), Step 1: design only, no code

Status: DRAFT for approval, 2026-10-12. Baseline code: git tag `baseline-pre-refactor`.
Anything marked **TBD** is not decided and needs a value from the owner.

## 0. Decisions this design follows (from the owner, not reopened)
1. Aligner stays STAR (+ featureCounts) with Salmon/tximport; no HISAT2 work.
2. One repo, two internal stages: QUANT and ANALYSIS. `nextflow run main.nf` with no `-entry` runs QUANT then ANALYSIS.
3. Significant gene = padj < 0.05 AND |apeglm-shrunk log2FC| >= 1. The unshrunk log2FC is always reported next to it.
4. Pathway background (ORA universe) = the genes DESeq2 actually tested, i.e. after the low-count pre-filter.
5. Sample sheet gets `count_source`, `baseline`, `batch`.
6. Two separate pre-filters stay, each a named parameter with today's value (section 6).
7. Pathway gene sets are read from local, versioned files whose name, version and md5 go into the manifest. No live web calls in the product path (this removes the KEGG-REST call and its MSigDB fallback in `ora.R`).
8. Validation-only scripts (`check_vs_paper.R`, `modules/check_vs_paper.nf`) leave the product path (section 8).
9. Group names and the baseline group come from the sample sheet, not from code.
10. Below 32 GB memory (or unreadable limit) QUANT stops with an error unless `--allow_salmon_only true`; the method actually used is written to the manifest and the report header.
11. KEGG stays on, behind `--run_kegg true|false`; licence for commercial use is unresolved and a note is printed in the report and docs.

## 1. Stages and entry points

| Entry | Runs | Needs | Produces |
|---|---|---|---|
| `-entry QUANT` | VERIFY_FASTQ, [SUBSAMPLE], FASTP, SALMON_QUANT, STAR_ALIGN, SAMTOOLS_*, FEATURECOUNTS, SAMPLE_QC_FLAGS, MERGE_*, SALMON_TXIMPORT, per-sample + cohort QC, MULTIQC (quant part), manifest | raw sample sheet + FASTQs + reference | the QUANT output contract (section 2) |
| `-entry ANALYSIS` | input validation, normalisation, DESeq2, shrinkage, ORA, GSEA, report | `--counts_matrix`, `--sample_sheet` (+ optional `--qc_table`, `--quant_manifest`) | DE tables, pathway tables, report, ANALYSIS manifest |
| default (no `-entry`) | QUANT, then ANALYSIS fed by QUANT's contract files | as QUANT | both |

- ANALYSIS never reads FASTQs, BAMs or `work/`; it reads only the contract files, so it can run on any matrix that satisfies the contract (for example the existing 19-sample matrix, or TCGA counts later).
- A selector `--analysis diagnostic` (default) is reserved so prognostic/predictive branches can be added later as sibling sub-workflows; only `diagnostic` exists now.
- Everything per sample stays per sample (no waiting for other samples); cohort steps still start after `collect()`.

## 2. QUANT -> ANALYSIS file contract

All files are tab-separated UTF-8, one header row, no quoting, written to `<outdir>/contract/` (the existing `01_`-`08_` folders and their contents are kept unchanged; the contract files are copies/derivations so the old layout and the regression comparison stay valid). Sample IDs are the SRA run accessions everywhere (existing naming rule); the label (T1_A ...) is an optional extra column, never a key.

### 2.1 `counts_matrix.tsv`
- Rows = genes (GENCODE v50 `gene_id`, versioned, first column `gene_id`), columns = sample IDs (SRR), values = **raw non-negative integers** from featureCounts (`-p --countReadPairs -B -C -Q 10`, strandedness per sample from Salmon).
- Column set = the samples that passed QC (not FAIL) and are in the sample sheet; gene order identical in all samples (already enforced by `merge_featurecounts.py`).
- No normalised, rounded or TPM values in this file.

### 2.2 tximport outputs (second route, optional)
- `counts_salmon_tximport.tsv`: tximport gene-level estimated counts (`countsFromAbundance = "no"`; non-integer, rounded only inside DESeq2 when built from a matrix) and `salmon_tpm.tsv`: gene TPM. Same row/column conventions as 2.1.
- `salmon_tximport.rds` is kept as the preferred ANALYSIS input for the Salmon route (it carries the length offset); `count_source = salmon_tximport` selects it.

### 2.3 `sample_sheet.tsv` (the one ANALYSIS reads)
| column | meaning |
|---|---|
| `sample_id` | SRR run accession, equals a column name in the matrix |
| `label` | optional human label (T1_A); never used as a key |
| `group` | condition name, any string; **names are not fixed in code** |
| `baseline` | `TRUE` for the reference group, `FALSE` otherwise; same value for all rows of a group; exactly one group is `TRUE` |
| `study` | source study / cohort |
| `batch` | batch factor for the design; if equal to `study` say so explicitly (for the current data `batch` = `study`) |
| `count_source` | `featurecounts` or `salmon_tximport`; how this sample's counts were produced |

- The contrast is every non-baseline group vs the baseline. Design = `~ batch + group` (today: `~ study + group` with `study` as batch, so results are unchanged when `batch = study`). More than two groups is **TBD**; v1 supports exactly two and stops with a message otherwise.
- The raw input sample sheet of QUANT (`sample,run,study,group`) gains the columns `batch`, `baseline`, `count_source` (QUANT fills `count_source` with the route actually used) and keeps working with the old 4-column format only if `--legacy_samplesheet true` (defaults: batch = study, baseline = group named `normal`). **TBD:** keep or drop this legacy switch.

### 2.4 `qc_table.tsv`
One row per sample: `sample_id`, `label`, `group`, `study`, `star_unique_pct`, `fc_assigned_pct`, `salmon_mapped_pct`, `fastp_q30_pct`, `fastp_pass_pct`, `status` (PASS|WARN|FAIL), **`reason`** (new: which metric crossed which limit, with the value and the limit; empty for PASS), plus cohort flags (`flag_outlier_mad`). Today's `sample_qc.tsv` has everything but `reason`. FAIL samples are excluded from the matrix and listed in the report with their reason. Limits come from the central thresholds file (Step 3).

### 2.5 `manifest.json`
Written by QUANT (`manifest_quant.json`) and by ANALYSIS (`manifest_analysis.json`, which embeds or references the former). Fields:
- `pipeline_version` (VERSION file / git tag), `git_commit` if available, `run_name`, start/end time, `entry`, `command_line`.
- `method_used`: `STAR+featureCounts` or `Salmon-only`; `allow_salmon_only`; cgroup memory limit read, and whether it was readable.
- `tool_versions`: Nextflow, Java, fastp, STAR, samtools, featureCounts, Salmon, R, DESeq2, tximport, apeglm, clusterProfiler, fgsea, msigdbr, etc.
- `parameters`: every resolved parameter value (all thresholds, section 6).
- `inputs`: md5 + size of sample sheet, counts matrix, QC table, every FASTQ (QUANT).
- `references`: md5 of genome FASTA, GTF v50, transcript FASTA, STAR index (`genomeParameters.txt` + SA/Genome md5 or size, **TBD** whole-index md5 is slow), Salmon index (`info.json`), tx2gene.
- `gene_sets`: for each database file: path, database name, release/version, source URL (informational only), date obtained, md5, number of sets.
- `sample_summary`: n per group, n per batch, n failed QC, baseline group name.
- `outputs`: md5 of the contract files.

## 3. ANALYSIS input checks (stop with a clear message, non-zero exit, nothing partial published)

Each check is a named function in one validation script with a one-line error message that says what was found and what is required.

1. **Baseline group missing / ambiguous**: no `baseline = TRUE` row, more than one baseline group, or `baseline` inconsistent within a group.
2. **Too few samples per group**: any group below `analysis_min_samples_per_group` (**TBD**; value to be agreed, no default assumed here) in the samples that pass QC. Message lists the group and counts.
3. **Counts content**: the integer check (non-integer values) applies **only when `count_source = featurecounts`**; tximport estimated counts (`salmon_tximport`) are legitimately non-integer and are not checked for integrality. For both sources: negative values, NA and duplicate gene IDs stop the run.
4. **Sample sheet and matrix mismatch**: any sample in one but not the other, duplicated IDs, or a column order that differs (reported as lists of the extra/missing IDs).
5. **`count_source` not recorded** or empty, an unknown value, or a mix of sources in one analysis.
6. **Design not estimable**: a batch x group cell empty (kept from the current `deseq2.R`), or batch perfectly confounded with group.
7. **Group/baseline names** read from the sheet; no hard-coded `normal`/`tumour`.
8. **Gene-set files**: any required gene-set file missing, or its md5 differing from the one in the manifest/config (no silent fallback, no download).
9. **Gene IDs**: matrix gene IDs must be found in the gene-id-to-symbol table; the fraction not found is reported and the run stops above a limit (**TBD**, parameter `analysis_max_unmapped_gene_pct`).

QUANT keeps its own existing guards (VERIFY_FASTQ gzip test, R1/R2 identical check, tximport label mapping) and adds the memory check below.

### 3.1 QUANT memory guard (main.nf)
- Read cgroup memory limit (v2, then v1). If unreadable, or < 32 GB: **error and stop**, unless `--allow_salmon_only true`, in which case STAR/featureCounts are skipped and `method_used = Salmon-only` is written to the manifest and printed in the report header. The 32 GB value becomes a named parameter (`min_mem_star_gb`).
- The `run_validation.sh` preflight (60 GB / 16 cores) is unchanged and stays outside the pipeline.

## 4. Statistics (ANALYSIS)
- Normalisation and model: DESeq2 size factors; `~ batch + group`; reference = baseline from the sheet.
- Filter before testing: gene needs >= `de_min_count` counts in >= `de_min_samples` samples (section 6). The set of genes that survive is saved as `tested_genes.tsv` and is the ORA background.
- LFC shrinkage: `lfcShrink(type = "apeglm")`. Output table columns: `log2FoldChange` (unshrunk), `log2FoldChange_shrunk`, `lfcSE`, `stat`, `pvalue`, `padj`.
- Significant: `padj < de_padj` and `abs(log2FoldChange_shrunk) >= de_lfc` (volcano classes, up/down tables, ORA input and concordance use this rule; today they use the unshrunk value).
- Normalised counts, VST, volcano, MA, top-50 heatmap as now (the 15 and 50 label/heatmap counts become parameters).

## 5. Pathways (ANALYSIS)
- Databases, each reported separately, BH correction **within** each database (never pooled): GO BP, GO MF, GO CC (simplified), KEGG, Reactome, optional MSigDB Hallmark.
- Each database gets: ORA (up and down separately; background = `tested_genes.tsv`) and GSEA on the ranked list. Ranking metric parameter `gsea_rank_metric` = `stat` (Wald) by default, alternative `lfc_shrunk` (**TBD:** confirm default).
- GO simplification: `clusterProfiler::simplify` with a similarity cutoff parameter (**TBD**, value to be agreed).
- Gene sets come from local versioned files in `ref/genesets/` (GMT, one file per database, e.g. `go_bp.<release>.gmt`, `reactome.<release>.gmt`, `kegg.<release>.gmt`, `msigdb_hallmark.<release>.gmt`), built once by a separate documented script and not by the pipeline run. The manifest records each file's version and md5. GO annotation versions (org.Hs.eg.db / GO.db release) are recorded as well. Where today's code uses `enrichKEGG` (web) or `msigdbr` (package data), the new code reads the GMT instead; the existing `msig_helper.R` is not used in the product path. **TBD:** which Reactome and KEGG release to freeze and where the KEGG GMT comes from (the MSigDB legacy KEGG set used today is already local package data but its licence status is also open).
- `--run_kegg false` skips KEGG ORA and GSEA entirely and records `kegg: skipped` in the manifest and report. When on, the report and `docs/` carry the note: "KEGG pathway data licence for commercial use is unresolved."
- Per-database methods text with tool and database versions is written by Step 4.

## 6. Parameters and the central thresholds file (Step 3 preview)
One file `conf/thresholds.config`, included from `nextflow.config`; no number in a script. Values here are today's and are not changed by the refactor, except where noted:

| Group | Parameter (proposed name) | Current value | Note |
|---|---|---|---|
| DE | `de_padj` | 0.05 | |
| DE | `de_lfc` | 1 | now applied to the **shrunk** value |
| DE pre-filter | `de_min_count`, `de_min_samples` | 10, 5 | name 1 of 2 pre-filters; test profile overrides to 2, 2 |
| Cohort-QC pre-filter | `cohort_min_count`, `cohort_min_samples` | 10, 3 (`min(3, n)`) | name 2 of 2; today hard-coded in `cohort_qc.R:12` |
| Cohort QC | `cohort_outlier_mad` | 3 | hard-coded `cohort_qc.R:28` |
| Sample QC | `qc_star_warn/fail`, `qc_fc_warn`, `qc_salmon_warn/fail` | 70/50, 60, 70/50 | in config today |
| fastp | `fastp_*` (6) | unchanged | |
| Memory guard | `min_mem_star_gb` | 32 | hard-coded `main.nf:87` |
| Samples | `analysis_min_samples_per_group` | **TBD** | |
| Gene IDs | `analysis_max_unmapped_gene_pct` | **TBD** | |
| Pathways | `path_min_gs`, `path_max_gs`, `path_q`, `gsea_nperm`, `gsea_seed` | 10, 500, 0.05, 10000, 42 | seed hard-coded today |
| Pathways | `go_simplify_cutoff` | **TBD** | |
| Switches | `run_kegg` | true | |

The resolved values are written to `manifest.json` (`parameters`).

## 7. What changes in existing files (for planning; nothing is edited in Step 1)
- `nextflow/main.nf`: split into `workflow QUANT`, `workflow ANALYSIS`, unnamed workflow = both; memory guard becomes an error; group/baseline from sheet (`['tumour','normal']` check removed).
- `deseq2.R`: baseline/group from sheet, `~ batch + group`, shrunk cutoff, `tested_genes.tsv`; plus a new validation script and a manifest writer.
- `ora.R`, `gsea.R`: local GMT, all databases, per-database BH, background from `tested_genes.tsv`; `msig_helper.R` and the web fallback leave the product path.
- `cohort_qc.R`: pre-filter and MAD from parameters. `sample_qc.py`: adds `reason`.
- `sample_map.nf` / sample-sheet parsing: new columns.
- Backups `.bak_pre_refactor` of every edited file; every change logged in `docs/run_log.md` and `docs/code_sources.md`.

## 8. Validation-only code leaves the product path
- Moved to `validation/` (new top-level folder, own entry `-entry VALIDATE_PAPER` or a plain script, **not** called by QUANT, ANALYSIS or the default run): `nextflow/bin/check_vs_paper.R`, `modules/check_vs_paper.nf`, and the Hosseini-specific gene/pathway lists inside it.
- `concordance.R` (featureCounts vs Salmon) compares two routes, not a paper; **TBD:** keep as an optional ANALYSIS output when both count sources exist, or move to `validation/` too (proposal: keep optional, off by default).
- Nothing from `validation/` feeds back into cutoffs (rule already in the script header).

## 9. Backward compatibility and what the regression test will compare (Step 5 preview)
- The default run keeps the `01_`-`08_` folders; new files are added (`contract/`, `07_pathway/<database>/`), none removed. `check_vs_paper_*` outputs leave the main `06_dge` folder (listed as an intended change).
- Intended differences vs the old 19-sample DESeq2 results: shrunk |LFC| cutoff (2,798 vs 3,213 genes on featureCounts, from the verification step), new gene sets/databases, local instead of web KEGG, ORA background = tested genes. Everything else must be identical; the changes are listed gene by gene.

## 10. Open items needing an owner decision (all marked TBD above)
1. `analysis_min_samples_per_group` value.
2. `gsea_rank_metric` default (Wald stat proposed).
3. `go_simplify_cutoff` value.
4. Reactome and KEGG release to freeze; KEGG GMT source.
5. Whether `analysis_max_unmapped_gene_pct` exists and its value.
6. Legacy 4-column sample sheet: keep switch or require new columns.
7. More than two groups: out of scope for v1 (proposal), confirm.
8. `concordance.R` placement.
9. Index md5 in the manifest: whole-index hash (slow) or size + parameters file.
10. Note: `results/10_independent_validation/scripts/run_validation.sh` sits under `results/` and is therefore not in the git baseline; say if it should be copied to `scripts/` and versioned.

## 11. Implementation notes (Steps 2-3, 2026-10-13; deviations and provisional choices, nothing here changes the owner's decisions)
- Provisional parameter values so the code runs (all in conf/thresholds.config, to be confirmed): `analysis_min_samples_per_group` = 3, `analysis_max_unmapped_gene_pct` = 5. `legacy_samplesheet` defaults to true (the existing diagnostic_20.csv has 4 columns).
- Salmon-only mode (`--allow_salmon_only true`): `contract/counts_matrix.tsv` is then the tximport counts TSV (non-integer), `count_source = salmon_tximport`; the integer check is not applied.
- Default run with the STAR route also runs the Salmon/tximport DESeq2 as a secondary cross-check (`SECONDARY_ROUTE`, validated with `--secondary true`, which accepts a sheet that records featurecounts) plus the concordance table, as before.
- ANALYSIS input checks: the QC table, when given, must not list FAIL samples that are still in the matrix; the sheet rows must be in matrix column order. The validated sheet and the check JSON are published to `00_analysis_input/` (`primary_*`, `secondary_*`).
- `-entry ANALYSIS` writes `RESULTS_INDEX_ANALYSIS.md` (it does not overwrite the QUANT `RESULTS_INDEX.md` in a shared outdir). The fastp QC copy in folders 06/07 comes from QUANT's MULTIQC step; an ANALYSIS-only run has no fastp input.
- Item 9 decision applied: manifest records md5 of `genomeParameters.txt` plus file sizes of the STAR index (no whole-index md5). FASTQ md5 in the QUANT manifest is on by default (`manifest_fastq_md5`).
- Left in nextflow/nextflow.config because they are resources/geometry, not thresholds: executor/process memory and cpus, `star_index_ram`, PAR coordinates, `wait_hours`; the 30 s FASTQ poll interval in modules/verify_fastq.nf is still a literal.
- ReactomePA is not installed in env/prefix/r_pathway; Step 4 reads Reactome from a local GMT anyway (design section 5).
