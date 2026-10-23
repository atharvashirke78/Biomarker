# Domain 1 of 4: Public end-to-end RNA-seq workflows vs ours

Review date: 2026-10-23. Reviewer: Claude agent (web tools: WebSearch, WebFetch, PubMed MCP).

## 0. Protocol (written BEFORE searching)

Question: which publicly available end-to-end RNA-seq workflows (FASTQ -> counts -> QC -> DE -> report) exist, and how do inputs, outputs, QC, reproducibility, licence and maintenance compare with ours (Nextflow DSL2, human PE bulk, fastp, STAR+featureCounts and Salmon/tximport, per-sample PASS/WARN/FAIL, cohort QC, DESeq2 ~batch+group + apeglm, ORA/GSEA, MultiQC, sha256 manifest, conda envs, no containers).

Sources and planned search strings (exact strings and counts are filled in section 1 as they are run):
1. GitHub (via WebSearch site:github.com): "RNA-seq pipeline nextflow DESeq2", "RNA-seq snakemake workflow DESeq2 differential expression", "bulk RNA-seq end-to-end pipeline FASTQ to differential expression report", "RNA-seq WDL pipeline STAR DESeq2".
2. nf-co.re pipeline list (nf-co.re/pipelines): browse for "rnaseq"/"differential"/"transcript".
3. Snakemake workflow catalog (snakemake.github.io/snakemake-workflow-catalog): search "rna-seq".
4. Galaxy Training Network / ToolShed: "reference-based RNA-seq", "RNA-seq workflow".
5. WorkflowHub: "RNA-seq".
6. Dockstore: "rnaseq".
7. bio.tools: "RNA-seq differential expression workflow".
8. PubMed / Europe PMC: "RNA-seq pipeline reproducible workflow differential expression", "Nextflow RNA-seq pipeline", "Snakemake RNA-seq".
9. bioRxiv: same keywords (note: bioRxiv MCP has no keyword search; web search used).

Inclusion: open-source licence; runs FASTQ (or counts) through to DE; documented; release or commit within the last 3 years (since 2023-10-23).
Exclusion: closed-source/commercial only; abandoned (no release/commit since before 2023-10-23); single-tool wrappers (e.g. only an aligner or only a DESeq2 wrapper with no upstream), non-RNA-seq (scRNA, ChIP) pipelines; duplicates.
Target: screen at least 40 candidates; log identified / screened / included / excluded with reasons (PRISMA-style).
Extraction per included workflow: engine, inputs, outputs, QC modules, DE methods, container/conda support, licence, last release date, URL.
Evidence marks: V = page opened and read; S = search snippet only. Failed pages listed explicitly.
Limitations declared in advance: web fetch summarises pages through a small model, so V facts are as reported by that summary; no workflows were executed; "last release" dates reflect page state on the review date.

## 1. Search log (what was actually run, 2026-10-23)

Access notes: WebFetch returns a small-model summary of a page, so "V" means the page (or raw file / API / feed) was opened and read, not that every sentence was verified word for word. The GitHub REST API unauthenticated limit (60/h) was exhausted after the first searches; afterwards raw.githubusercontent.com files and github.com `releases.atom` / `commits.atom` feeds were used for release and commit dates (counted as V, the feed was opened). No workflow was executed.

| Source | Exact string / URL used | Result |
|---|---|---|
| nf-co.re | `https://nf-co.re/pipelines.json` (the HTML list page is JS-rendered and returned no names; failed, see section 6). Filter: name/description/topics contain rna, transcript, expression or differential | 159 pipelines in catalogue; 45 keyword matches |
| GitHub (repo search API) | `rna-seq nextflow deseq2` (51 hits, top 30 read); `rna-seq snakemake deseq2` (53 hits, top 30 read); `rnaseq pipeline differential expression fastq report` (1 hit) | 61 records read, 60 unique |
| GitHub via WebSearch | `RNA-seq pipeline nextflow DESeq2 differential expression github`; `RNA-seq snakemake workflow DESeq2 differential expression github` | pointers to qbic-pipelines/rnadeseq, BioinformaticsSolutionCenter fork, niekwit, WorkflowHub 1963 |
| Snakemake workflow catalog | `https://snakemake.github.io/snakemake-workflow-catalog/docs/all_standardized_workflows.html` (the `?topics=rna-seq` page is JS-rendered; failed). Filter on names containing rna, expression, deseq, transcript | 45 name matches |
| WorkflowHub | `https://workflowhub.eu/workflows?filter[query]=RNA-seq&format=json` | 119 records |
| Dockstore | `https://dockstore.org/api/workflows/published?search=rnaseq&limit=100` (first attempts without limit or at `/entries/search` failed) plus WebSearch `Dockstore RNA-seq workflow WDL STAR RSEM differential expression` | 100 records read (total unknown, limit is 100); only nf-core/rnaseq and a HuBMAP Salmon WDL are bulk RNA-seq |
| Galaxy ToolShed | `https://toolshed.g2.bx.psu.edu/api/repositories?q=deseq2` | 15 hits, all tool wrappers (single-tool, excluded as a class) |
| Galaxy Training / GTN | WebSearch `Galaxy Training Network reference-based RNA-seq data analysis tutorial DESeq2 workflow`; opened `training.galaxyproject.org/.../ref-based/tutorial.html` | 1 tutorial + workflows (WorkflowHub 1711, 1713, 1715) |
| bio.tools | `https://bio.tools/api/tool/?q=RNA-seq differential expression workflow&format=json` | 96 hits; first 35 read (rest not read, page truncated) |
| PubMed | `RNA-seq pipeline reproducible workflow differential expression Nextflow OR Snakemake` (312 hits, 20 titles read); `Nextflow RNA-seq pipeline` (48 hits, 20 read); `Snakemake RNA-seq workflow differential expression` (15 hits, titles read for the first ~21 of the 33 unique PMIDs in the second and third queries) | 1 relevant new item (ARACRA, 2026) plus RNA-SeqEZPZ |
| Europe PMC | NOT searched (no tool; PubMed MCP used instead) | - |
| bioRxiv | The bioRxiv MCP has no keyword search. WebSearch `bioRxiv end-to-end bulk RNA-seq pipeline FASTQ differential expression automated report Nextflow` | no relevant preprint found; not a systematic bioRxiv search |
| Named by knowledge or web snippets, then checked via feeds/READMEs | ENCODE rna-seq-pipeline, bcbio, snakePipes, RENEE, ZARP, GEMmaker, RNAflow, Curare, ARMOR, RASflow, Sequana, qbic rnadeseq, nextflow-io/rnaseq-nf, ARACRA, HULK | see table |

## 2. PRISMA-style flow (counts are unique candidate workflows after merging duplicates across sources)

```
Records read from sources (before de-duplication)
  nf-core 45 | GitHub search 61 | Snakemake catalog 45 | WorkflowHub 119 (26 RNA-bulk-relevant after title screen)
  | bio.tools 35 read | Dockstore 100 read (2 relevant) | ToolShed 15 | GTN 1 | PubMed ~60 titles read | named 15
Unique candidate workflows screened at record level ............ 202
   nf-core 45 + GitHub (new) 59 + Snakemake catalog (new) 41 + WorkflowHub (new) 26 + other named/bio.tools/Dockstore 31
Excluded ............................................................. 181
Included ............................................................. 21
```

Exclusion reasons (181, each candidate counted once; sums checked):
- out of scope (not bulk short-read gene-level RNA-seq to DE: single-cell, long-read, small RNA, fusion, splicing, ChIP, variant calling, spatial, bacterial, epitranscriptome, teaching, non-pipeline): 72 (nf-core 43, GitHub 3, Snakemake catalog 25, WorkflowHub 1)
- no licence file in the repository record: 30 (all GitHub)
- listed but not opened (description line or title only, single-author or single-dataset repos, or no code found): 41 (GitHub 13, Snakemake catalog 9, WorkflowHub 12, other 7: ARACRA, DEGgo, RseqFlow, SPARTA, aTAP, Myrna, RNASSIST). These are NOT judged bad, only unassessed.
- stops before DE (FASTQ to counts only): 16 (Snakemake catalog 5, WorkflowHub 4 incl. ZARP and IWC rnaseq-pe/sr, other 7 incl. RENEE, GEMmaker, Sequana, nextflow-io/rnaseq-nf)
- stale or discontinued: 7 (GitHub 4 with licence, last push 2019-2021; RASflow 2021; ENCODE rna-seq-pipeline last commit 2022-03-31; bcbio-nextgen discontinuation notice 2024-08-16)
- duplicate, fork or sibling of an included entry: 5 (one further duplicate GitHub record was removed before counting)
- R package, web app or tool wrapper (DE only, no FASTQ stage): 5 (VISTA, KnowSeq, compcodeR, CoRMAP, Galaxy ToolShed DESeq2 wrappers)
- closed source (Partek Flow, CLC, DRAGEN, BaseSpace): 4, named from prior knowledge, not searched or opened
- repository not resolvable (kundajelab/rna-seq-pipeline): 1
Per-source figures (45 nf-core, 59 GitHub new, 41 Snakemake new, 26 WorkflowHub new, 31 other = 202) were counted by hand from the lists read; the per-source bullets in section 3 are the authoritative record.

## 3. Screening table (summary; I = included, X = excluded)

nf-core catalogue (45 keyword matches): I rnaseq (conditional: DE is not inside it, see section 4, row 1), I differentialabundance; X 43 others, none is a bulk gene-level FASTQ-to-DE workflow: alleleexpression, ampliseq, cageseq (2021), chipseq, circrna, clipseq (2021), daamicrobiome, dartseq, denovotranscript, drop, dualrnaseq (2021), evexplorer, hicar, hlatyping, isoseq, lncpipe, marsseq, metatdenovo, methylarray, molkart, ncrnannotator, nanostring, nascent, oncoanalyser, rnadnavar, rnafusion, rnasplice (splicing), rnastructurome, rnavar, scdownstream, scflow, scnanoseq, scrnaseq, slamseq (2020), smartseq2 (2021), smrnaseq, sopa, spatialaxe, spatialvi, stableexpression, tfactivity (V: pipelines.json).

GitHub search (60 unique): I snakemake-workflows/rna-seq-star-deseq2, tucca-cellag/tucca-rna-seq, niekwit/rna-seq-salmon-deseq2, schlesnerlab/multiconDEnrich, challagandla/oracle-rna-seq, mreyer93/bulk-rnaseq-pipeline, Millimono/OmicsFlow, gyagu98/rna-seq-nextflow-pipeline, oxo-flow-community/oxo-flow-rnaseq-star-deseq2 (nf-core/differentialabundance counted above). X (50): 30 without any licence in the API record (e.g. ashlesha874/nextflow-rnaseq-pipeline, DimAnge/bulk-rnaseq-lite-nf, barbavegeta, shruthisundar793, kGorze, pandeyravi15, and 24 more), 4 stale with licence (MLKaufman 2021, TranslationalBioinformaticsUnit 2019, arabidopsisca 2020, moritzschaefer 2019), 3 out of scope (challagandla/oracle-mirna-seq, hmkim migration guide, josh45 re-analysis), 13 MIT/GPL single-project or thin-wrapper repos not opened (Muyi-bioinfo, alexturcoo, VladNikitin, Cslaer2025, Aabida1, michalbukowski, 7-jayanth2003, xHaMMaDy, niekwit/rna-seq-star-deseq2 sibling, GenomeEditingPlatformHZAU, umranyaman, Kin-zala, tdayris/fair_rnaseq_deseq2 whose README returned 404), 1 duplicate record. (S for the 13: description line only.)

Snakemake catalog (45, 4 duplicates of the above, 41 new): I snakemake-workflows/rna-seq-kallisto-sleuth; X 40 (no DE: Bulk-RNA-Snakeline V, epigen/rnaseq_pipeline V, baerlachlan counts, tdayris salmon_quant, mpitz123; out of scope: long-read, bacterial, small RNA, single-cell, variant calling, epitranscriptome, 25 in all; not opened: ArsenaultResearch, FeelLiao/rna-seq-std, bakerwm, gynecoloji, jma1991, loganylchen/RNA-seq-analysis-pipeline, osvaldoreisss, rochevin, vjcitn; fork: plycrsk).

WorkflowHub (26 new relevant titles): I 1201 IWC rnaseq-de, 1713 GTN DEG, 1814 RNA-SeqEZPZ-NF, 1963 RNAseq-quantif-pipeline; X 22 (400/401 IWC rnaseq-sr/pe: FASTQ to counts, no DE, V for 401; 2278 oxo-flow port of nf-core/rnaseq: no DE, V; 447 ZARP no DE; 451 RASflow stale; 1813 duplicate UI of 1814; 2300 sibling of 2290; 1711/1715 GTN sub-workflows; 1688 isoform scope; 58, 524, 152, 411, 412, 685, 689, 1966, 1692, 1699, 1703, 338 titles only, S).

Other named (31): I snakePipes mRNAseq, RNAflow, Curare, ARMOR, qbic-pipelines/rnadeseq; X 26: ENCODE-DCC/rna-seq-pipeline (last commit 2022-03-31, last release 2021), kundajelab/rna-seq-pipeline (no repo resolved), bcbio-nextgen (README: "notice of discontinuation ... 08-16-2024", last commit 2024-08-24), RENEE (stops at counts), GEMmaker (last release 2023-02-03, matrix only), Sequana rnaseq (counts only), nextflow-io/rnaseq-nf (demo, no DE), ARACRA (PubMed 2026; no code repository found, so cannot be assessed), HULK (S, no DE), HuBMAP Salmon WDL, TOPMed CWL, DEGgo, RseqFlow, SPARTA, aTAP, Myrna, RNASSIST (bio.tools snippets only), VISTA, KnowSeq, compcodeR, CoRMAP (R packages), ToolShed DESeq2 wrappers, Partek Flow, CLC, DRAGEN, BaseSpace (closed).

## 4. Extraction for the 21 included workflows

V = page, raw file, feed or JSON opened; S = snippet or title only. Dates are the latest GitHub release (feed) unless noted; "commit" = latest commit in the commits feed. "DE" = differential expression.

| # | Workflow | Engine | Inputs | Outputs and QC modules | DE methods | Env support | Licence | Last release / commit | URL |
|---|---|---|---|---|---|---|---|---|---|
| 1 | nf-core/rnaseq 3.27.0 (conditional: no DE inside, pairs with #2) | Nextflow | CSV: sample, fastq_1, fastq_2, strandedness (`auto` infers with Salmon); optional UMI, rRNA removal (V) | Counts and tximport matrices (SummarizedExperiment), bigWig, MultiQC. QC: fq lint (only hard stop), FastQC, RSeQC (7 scripts), Qualimap, dupRadar, Preseq, featureCounts biotypes, Kraken2/Bracken/Sylph, strandedness and low-mapping flags in MultiQC (flag only, no sample exclusion), experimental RustQC (V) | None. DESeq2 only for PCA/sample-distance QC; docs send users to #2 (V) | docker, singularity, apptainer, podman, wave; conda "last resort" (V) | MIT (V) | 2026-09-23 / commit 2026-10-08 (V); 1394 stars (V pipelines.json) | https://nf-co.re/rnaseq/3.27.0/docs/output/ , /usage/ |
| 2 | nf-core/differentialabundance 2.0.0 | Nextflow | Samplesheet (batch/blocking columns), abundance matrix, contrasts CSV/YAML, GTF (V) | DE tables, Quarto report, ShinyNgs app (V). No per-sample QC gate described | DESeq2 (default for RNA-seq), limma, DREAM (random effects), propd; interactions; GSEA, gProfiler2, decoupler (V). LFC shrinkage not stated on the page read | Container profiles; conda last resort (V) | MIT (V) | 2026-06-23 / commit 2026-10-08 (V); 95 stars | https://nf-co.re/differentialabundance/2.0.0/docs/usage/ |
| 3 | snakemake-workflows/rna-seq-star-deseq2 v3.1.2 | Snakemake >= 8 | samples.tsv + units.tsv (lanes), local FASTQ or SRA accession, Ensembl species/release auto-download, strandedness per unit (V) | fastp (since v3.0.0; HTML in Snakemake report), STAR GeneCounts, PCA per variable (V). No MultiQC or gating seen in the files read | DESeq2, model `~batch_effects + variables_of_interest` (default includes interactions), contrasts (V config). r-ashr is in the DESeq2 env, so ashr shrinkage is probable (S). No enrichment | conda pinned per rule (V); containers via Snakemake (S) | MIT (V) | 2026-10-07 (CHANGELOG, V); 373 stars; Zenodo DOI; CI (V) | https://github.com/snakemake-workflows/rna-seq-star-deseq2 |
| 4 | snakemake-workflows/rna-seq-kallisto-sleuth v4.1.1 | Snakemake | per catalog (S) | kallisto; Snakemake report (S) | sleuth (V description) | not stated | MIT (V) | 2026-10-15 (feed). The page also showed an "Archived" label that conflicts with a 2026 release; unresolved | https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth |
| 5 | niekwit/rna-seq-salmon-deseq2 v0.9.0 | Snakemake >= 8.25 | PE and SE FASTQ auto-detected, samples.csv (genotype, treatment, reference, batch), GENCODE human/mouse auto-download (V) | FastQC, trimming, Salmon, volcano plots (V) | DESeq2 at gene and transcript level, default `~batch + comb`, custom design (V). Shrinkage not stated | Apptainer container with conda envs baked in, or conda (V) | MIT (V) | 2026-09-22 (feed); Zenodo DOI; CI | https://github.com/niekwit/rna-seq-salmon-deseq2 |
| 6 | tucca-cellag/tucca-rna-seq v1.0.1 | Snakemake >= 8.27 | raw FASTQ, metadata, JSON-schema-checked config, contrasts, optional GMT (V) | FastQC, Qualimap, Salmon, MultiQC, GeneTonic and pcaExplorer Shiny apps (V) | DESeq2 (rlog/vst; no shrinkage named); clusterProfiler, SPIA, MSigDB enrichment (V) | conda, singularity, apptainer badges (V) | MIT in GitHub listing (V); the docs page itself states none | 2026-03-02 (feed); docs say it is being reworked as a downstream toolkit for nf-core/rnaseq | https://github.com/tucca-cellag/tucca-rna-seq , https://tucca-compbio.tufts.edu/tucca-rna-seq/introduction |
| 7 | schlesnerlab/multiconDEnrich v0.2.0 | Snakemake | counts (alignment removed from #3), multi-group config | HTML per comparison (V) | DESeq2, enrichment, optional CARNIVAL (needs CPLEX solver) (V) | conda | MIT (S, search API) | 2024-10-31 (feed) | https://github.com/schlesnerlab/multicondition-deseq2-enrichment |
| 8 | challagandla/oracle-rna-seq (closest analogue to ours) | Snakemake | PE FASTQ, samples.tsv, human/mouse/rat; strict validation of IDs, pairs, replicates, design rank (V) | FastQC, fastp, optional FastQ Screen, STAR two-pass, Picard, RSeQC, Qualimap, preseq, strandedness inference that stops the run if ambiguous, featureCounts (-B -C), Salmon, StringTie, MultiQC, resolved config saved with results (V) | DESeq2 with shrinkage that must succeed (type not stated), GO/KEGG ORA on tested-gene universe, GSEA ranked by Wald statistic (V) | 17 isolated conda envs, checksum-verified Miniforge; no containers (V) | MIT own code, third-party licence notice (V) | no release (feed empty); pushed 2026-07-15 (API) | https://github.com/challagandla/oracle-rna-seq |
| 9 | mreyer93/bulk-rnaseq-pipeline | Snakemake | sample sheet with flexible column names, line-numbered validator; SE or PE; lanes merged (V) | fastp, FastQC, MultiQC, Salmon (decoy-aware) or STAR+Salmon, tximport, HTML/PDF report (V) | DESeq2 `batch + condition`, contrasts (V) | conda (`--use-conda`) (V) | Apache-2.0 (V) | no release; pushed 2026-10-02 (API) | https://github.com/mreyer93/bulk-rnaseq-pipeline |
| 10 | Millimono/OmicsFlow v1.0.2 | Nextflow | FASTQ; STAR index | FastQC, STAR, Salmon, volcano, heatmap, PCA; validated on one yeast test set (V) | DESeq2 (V) | Docker image 4.6 GB (V) | MIT (V) | 2026-06-23 (feed); DOI | https://github.com/Millimono/OmicsFlow |
| 11 | gyagu98/rna-seq-nextflow-pipeline | Nextflow DSL2 | samplesheet.csv | FastQC, Salmon, MultiQC, CI with stub run (V) | DESeq2 (V) | Docker, Singularity, conda (V) | MIT (V) | no release; pushed 2026-10-22 (API) | https://github.com/gyagu98/rna-seq-nextflow-pipeline |
| 12 | oxo-flow-community/oxo-flow-rnaseq-star-deseq2 (WorkflowHub 2290) | oxo-flow (typed "WDL" on WorkflowHub) | units sheet with SRA option, per-unit fastp overrides (V) | fastp, STAR, RSeQC, MultiQC, biomaRt annotation, PCA, MA plots (V) | DESeq2 per contrast with ashr shrinkage (V) | every tool pinned to a conda version (V) | "Other (Open)" on WorkflowHub; Apache-2.0 in GitHub listing | frozen 2026-09-05, single submitter, new and unproven | https://workflowhub.eu/workflows/2290 |
| 13 | IWC rnaseq-de v0.12 | Galaxy | two collections of count files + GTF, padj and LFC thresholds (V) | Annotated DESeq2 table, volcano, two heatmaps, normalised counts (V) | DESeq2; two-condition designs only (V) | Galaxy tool dependencies | MIT (V) | 2026-09-07 (V WorkflowHub; feed) | https://workflowhub.eu/workflows/1201 , https://github.com/iwc-workflows/rnaseq-de |
| 14 | Galaxy Training Network "Reference-based RNA-seq" + workflows | Galaxy | FASTQ (SE+PE), Drosophila example (V) | Falco/FastQC, MultiQC, Cutadapt, STAR, featureCounts, RSeQC, Picard (V) | DESeq2, goseq, KEGG (V from step list; later sections not read) | Galaxy servers | CC-BY 4.0 tutorial, MIT framework (V) | tutorial modified 2026-09-23 (V); WorkflowHub 1713 v8 2024-07-05 (S) | https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html |
| 15 | RNA-SeqEZPZ-NF v1.0.1 | Nextflow | PE FASTQ, genome FASTA, GTF via a Shiny form (V) | Trimming, QC tab, SARTools-based statistical report, pathway analysis, UpSet/volcano (V) | DESeq2 (V) | Singularity image; conda not mentioned (V) | LICENSE file present, type not read | 2025-07-22 (feed); 0 stars; paper in PubMed | https://github.com/yzhang18/RNA-SeqEZPZ-NF |
| 16 | RNAseq-quantif-pipeline (WorkflowHub 1963) | Snakemake | PE only | prinseq/cutadapt, STAR, HTSeq (V) | DESeq2, clusterProfiler GO/KEGG (V) | not stated | GPL-3.0 (V) | 2025-09-25 (V) | https://workflowhub.eu/workflows/1963 |
| 17 | snakePipes mRNA-seq v3.5.0 | Snakemake | FASTQ dir (SE/PE) or BAM, sample sheet with condition and blocking columns (V) | fastp/cutadapt/TrimGalore, STAR or HISAT2, featureCounts, deepTools QC, MultiQC, DESeq2 HTML report, allelic modes (V) | DESeq2 (`--formula`, `--LRT`), salmon + sleuth/wasabi, rMATS (V) | conda by default; containers not described (V) | not verified (page silent, no LICENSE found at probed paths) | 2026-10-21 (feed) | https://snakepipes.readthedocs.io/en/latest/content/workflows/mRNAseq.html |
| 18 | hoelzer-lab/RNAflow 1.6.1 | Nextflow >= 20.10 | reads CSV (condition, source, strandedness), auto-download genomes (V) | fastp, SortMeRNA, HISAT2, featureCounts, TPM filter, MultiQC, MA/PCA/heatmap/volcano (V) | DESeq2 all pairwise by default; WebGestalt and piano GSEA (V) | conda, docker, singularity (V) | GPL-3.0 (V) | 2025-07-07 (feed); commit 2025-08-06 | https://github.com/hoelzer-lab/rnaflow |
| 19 | Curare 0.6.0 | Snakemake + conda | samples.tsv + pipeline.yml via wizard (V) | Trim Galore/fastp, FastQC, MultiQC, many mappers, featureCounts, coverage, HTML and Excel (V) | DESeq2 (V) | conda only; `--no-conda` option (V) | GPL-3.0 (V) | 2024-02-27 (feed) | https://github.com/pblumenkamp/Curare |
| 20 | csoneson/ARMOR v1.5.11 | Snakemake | FASTQ, config.yaml | QC, optional STAR, Salmon, DRIMSeq, iSEE app; details thin (V) | DESeq2/edgeR (page names them only as swappable) (V) | conda (V) | MIT (V) | 2026-06-14 (feed); 167 stars | https://github.com/csoneson/ARMOR |
| 21 | qbic-pipelines/rnadeseq 2.6 | Nextflow (nf-core template, >= 19.04) | counts + metadata (downstream only) | DESeq2 reports (S for details) | DESeq2, pathway analysis (V description) | Docker, bioconda (V) | MIT (V) | 2025-07-03 (feed); 32 stars (S) | https://github.com/qbic-pipelines/rnadeseq |

## 5. Comparison with ours

Ours (from `README.md`, `conf/thresholds.config`, user brief): Nextflow DSL2, human paired-end bulk RNA-seq only, fastp, STAR + featureCounts and Salmon/tximport, per-sample PASS/WARN/FAIL (FAIL excluded from cohort steps; thresholds on STAR, featureCounts and Salmon mapping %), cohort QC (PCA, correlation, library size, MAD flags), DESeq2 `~batch + group` with apeglm, ORA and GSEA (fgsea, seed fixed), MultiQC, sha256 manifest, `params.json`, `tool_versions.tsv`, conda envs, no containers, no LICENSE file yet. "No evidence" below means not seen in the pages I read, not proof of absence.

| Feature | Ours | nf-core rnaseq + diffabund | snakemake star-deseq2 | ORACLE | mreyer93 | RNAflow | snakePipes | tucca |
|---|---|---|---|---|---|---|---|---|
| Engine | Nextflow DSL2 | Nextflow | Snakemake | Snakemake | Snakemake | Nextflow | Snakemake | Snakemake |
| Layout and species | human PE only | SE/PE, any species | SE/PE (SRA PE), Ensembl species | PE, human/mouse/rat | SE/PE | SE/PE, auto-download | SE/PE | not stated |
| Lanes per sample | not supported | merged | units sheet | not stated | merged | not stated | not stated | not stated |
| Trimming | fastp | Trim Galore or fastp | fastp | fastp | fastp | fastp | choice of 3 | not stated |
| Quantification | STAR+featureCounts and Salmon/tximport | STAR+Salmon (default), RSEM, kallisto | STAR GeneCounts | featureCounts and Salmon | Salmon or STAR+Salmon | HISAT2+featureCounts | featureCounts or Salmon | Salmon |
| Per-sample gate with exclusion | yes (PASS/WARN/FAIL with reason) | no: flags only, fq lint hard stop | no evidence | strand ambiguity stops run; validation | validation only | no evidence | no evidence | no evidence |
| Depth of QC metrics | mapping %, cohort PCA/MAD; RSeQC/dupRadar/Qualimap not found in our modules by grep (not exhaustively checked) | RSeQC, Qualimap, dupRadar, Preseq, Kraken2, strandedness check | fastp only | RSeQC, Picard, Qualimap, preseq, FastQ Screen | FastQC, fastp | fastp, SortMeRNA | deepTools, fastp | FastQC, Qualimap |
| Cohort QC | PCA, correlation, MAD flags | DESeq2 PCA, sample distance | PCA per variable | not stated | not stated | PCA, heatmap | PCA, correlation | pcaExplorer |
| Design | `~batch + group` | formula, interactions, random effects, blocking | `~batch + variables` with interactions | design + contrasts | `batch + condition` | pairs + source | `--formula`, LRT | designs + contrasts |
| LFC shrinkage | apeglm | not stated | ashr in env (S) | yes, type not stated | not stated | not stated | not stated | not named |
| Enrichment | ORA + GSEA | gProfiler2, GSEA, decoupler | none | GO/KEGG ORA, GSEA | none | WebGestalt, piano | rMATS | clusterProfiler, SPIA |
| Report | MultiQC + HTML/MD report | MultiQC, Quarto, Shiny | Snakemake report | MultiQC | HTML/PDF | MultiQC + HTML tables | DESeq2 HTML, MultiQC | MultiQC + Shiny |
| Provenance | sha256 manifest of inputs, configs, lock, outputs; run_info, params, tool_versions | software_versions, params, pinned releases | Snakemake report, Zenodo DOI | resolved config saved | per-rule logs | Nextflow report | log, DAG | JSON-schema config |
| Environments | conda only | containers, conda last resort | conda (+ containers via Snakemake) | conda only | conda | conda/docker/singularity | conda | conda, singularity |
| Licence | none yet | MIT | MIT | MIT | Apache-2.0 | GPL-3.0 | not verified | MIT |
| Public release and DOI | none (no git repo yet) | releases, DOIs | v3.1.2, Zenodo | none | none | 1.6.1 | v3.5.0 | v1.0.1, Zenodo |

### 5.1 Gaps we have
1. No licence file and no tagged public release or DOI; every mature peer (#1-3, #5, #6, #18, #20) has both. This blocks reuse and citation.
2. No container option. nf-core, RNAflow, niekwit, tucca and RENEE-style pipelines offer Docker/Apptainer; conda alone drifts as channels change (nf-core itself calls conda a last resort).
3. Narrow QC depth: no RSeQC (infer_experiment, read distribution), dupRadar, Preseq, rRNA or contamination screening seen in our modules; nf-core and ORACLE have them. Strandedness inference (nf-core `auto`, ORACLE confident-call rule) is not a feature I verified in ours.
4. Input scope: no lane merging, single-end, other species, or SRA input (nf-core, #3, #9).
5. Design flexibility: only `~batch + group`; peers offer contrasts files, interactions, blocking, random effects (nf-core diffabund, snakePipes LRT). Our README lists optional `patient_id`; whether a paired `~patient + group` design is used was not checked and matters for tumour-normal colon data.
6. Alternative DE engines and transcript-level DE (limma/DREAM in #2, sleuth in #4/#17, transcript level in #5).
7. No interactive report apps (Shiny in #2, #6, #15); low priority for a locked-model pipeline.

### 5.2 Strengths we have
1. Explicit per-sample PASS/WARN/FAIL with a named reason and automatic exclusion of FAIL samples from cohort steps. In the pages read, nf-core only flags (MultiQC highlight, no exclusion), ORACLE stops on ambiguity, others show nothing comparable.
2. All thresholds in one file (`conf/thresholds.config`) written into the run manifest; few peers document this.
3. sha256 manifest covering inputs, configs, the locked model and outputs, plus `params.json` and `tool_versions.tsv`. No peer page I read describes output checksums (Snakemake peers rely on reports and DOIs; nf-core on pinned releases and `software_versions`).
4. apeglm shrinkage stated explicitly; peers state ashr (#12), none (#2 page read), or do not say.
5. Two quantification routes in one run (featureCounts and Salmon/tximport) with GSEA seed fixed; locked-model apply mode with hash check, which none of the 21 offers.
6. Honest "NOT TESTED" labelling and a stub test suite; most small peers (#9-#12) have only CI toy runs.
7. Single-command input scan with exit before compute (comparable to #8 and #9 validators).

### 5.3 Ranked adoptions (benefit over effort)
1. Add a LICENSE file, git history, a tagged release and a Zenodo DOI (hours; unblocks everything else; owner decision per README).
2. Optional container profile (Docker/Apptainer from BioContainers pins) next to conda, keeping conda as default; at minimum export explicit conda lock files (hours to days; copies nf-core profiles and #5).
3. Strandedness inference with a WARN/FAIL (Salmon `auto` as in nf-core, or RSeQC infer_experiment) plus read-distribution and duplication metrics as WARN-only (days; closes the biggest QC-depth gap without changing verdict rules).
4. Check and, if absent, add paired design `~patient_id + group` when pairs exist, and a contrasts/blocking option modelled on nf-core/differentialabundance (days; high scientific value for tumour-normal).
5. Cross-check our counts and DE on the 19 development samples against nf-core/rnaseq 3.27.0 + differentialabundance 2.0.0 (one run each; yields an external concordance figure; do not adopt their code).
6. Lane merging via a units sheet (#3, #1) (days).
7. Flag-only contamination screen (Kraken2 or FastQ Screen, as in #1, #8) (days).
8. Do not adopt: switching engine to nf-core wholesale (loses locked model, gating and manifest), interactive Shiny, long-read or splicing modules.

## 6. Pages that failed or were only partly read

- https://nf-co.re/pipelines : rendered text had no pipeline list; replaced by https://nf-co.re/pipelines.json (opened).
- https://snakemake.github.io/snakemake-workflow-catalog/?usage=all&topics=rna-seq : JS-rendered, no rows; replaced by `docs/all_standardized_workflows.html` (names only, no stars/dates, which I took from GitHub feeds).
- GitHub REST API: rate-limited (60/h) after search; repo metadata for #1-#21 came from feeds and raw files instead. Licence "S, search API" entries (#7, #12) came from the search listing only.
- https://tucca-cellag.github.io/tucca-rna-seq/introduction : 301 redirect, followed manually to tucca-compbio.tufts.edu (opened).
- https://snakepipes.readthedocs.io/.../mRNA-seq.html : 404, alternative `mRNAseq.html` opened. snakePipes licence not found (probed LICENSE, LICENSE.md, LICENSE.txt on master).
- https://ccbr.github.io/RENEE/ : redirect-only; `/latest/` opened (RENEE then excluded, no DE).
- https://dockstore.org/api/workflows/entries/search?query=rnaseq : 404; the `published?search=rnaseq&limit=100` endpoint worked but returns at most 100 records and the full result set was not paged.
- https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html : 404; README raw file used for the discontinuation notice.
- tdayris/fair_rnaseq_deseq2 README : 404. kundajelab/rna-seq-pipeline : no feed (repo moved or removed).
- Truncated by the fetch tool: bio.tools list (35 of 96 read), WorkflowHub pages 2290/2300/1201/1814/2278/1963 (last 5-20 KB not read), GTN tutorial (first 100 of 197 KB read, so DE and enrichment sections were not read beyond the step list).
- Not done: Europe PMC search, Galaxy ToolShed workflow search beyond `deseq2`, a systematic bioRxiv search, the full Snakemake catalog "all other workflows" list, GitHub code search, Bioconda/BioContainers cross-check, any execution of any workflow.
- ARACRA (PubMed 42620692, Nextflow + Streamlit, FASTQ to transcriptomic point of departure, published 2026-08-19) has no code repository I could find, so it could not be assessed or included.
- Conflicts left unresolved: WorkflowHub 401 shows "last updated 2023-01-17" next to version v1.6 dated 2026-10-15; snakemake-workflows/rna-seq-kallisto-sleuth page shows "Archived" next to a 2026-10-15 release feed entry; star counts and licence strings come from the moment of reading and from a model summary of the page.
- Closed-source tools (Partek Flow, CLC, DRAGEN, BaseSpace) were named from prior knowledge and not searched.

## 7. Bottom line

Public end-to-end options that go FASTQ to DE with reports exist, but only two are widely used and actively maintained at scale: nf-core/rnaseq paired with nf-core/differentialabundance (Nextflow, MIT, containers) and snakemake-workflows/rna-seq-star-deseq2 (Snakemake, MIT, 373 stars). Several 2025-2026 single-author Snakemake/Nextflow repos (#8, #9, #10, #11, #12) have the same shape as ours; the closest, ORACLE (#8), matches our scope but has no release and no per-sample gate. We are ahead on gating, threshold transparency and output provenance, and behind on licensing and release hygiene, container support, QC depth and design flexibility.
