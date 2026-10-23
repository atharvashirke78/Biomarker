# Review 03 of 4: RNA-seq QC checks and reproducibility practices

Date of search: 2026-10-23. Marking: **V** = verified by me this session against the primary page, file or API response. **S** = secondary, recalled, or inferred; not confirmed this session. "Threshold: none" means I looked in the cited source and it states no number.

## 1. Protocol (written before searching)

Question: which RNA-seq QC checks and reproducibility practices are recommended or standard (strandedness, rRNA, 3' bias/TIN, duplication, gene-body coverage, contamination, sex check, sample swap/identity, replicate correlation, library complexity; containers, version pinning, BLAS/thread determinism, Salmon non-determinism)? Which tools, with which thresholds that a source states, and under which licences?

Our current state: per-sample QC from STAR unique %, featureCounts assigned %, Salmon mapped %, fastp Q30; cohort PCA and MAD flags; flag-only sex check (XIST/Y genes); strandedness from Salmon library type only; conda envs; sha256 manifest.

Planned sources and exact strings (all run on 2026-10-23):

| Source | Access method | Exact search strings |
|---|---|---|
| PubMed (MCP tool, natural language) | `search_articles` | (1) `RNA-seq quality control metrics reproducibility RSeQC OR Qualimap OR dupRadar`; (2) `RNA-seq sample swap identity check genotype concordance OR "sex check" quality control`; (3) `transcript integrity number RNA-seq degradation 3' bias quality`; (4) `Salmon quantification non-deterministic OR reproducibility RNA-seq pipeline bioinformatics software versions containers`; (5) `contamination detection RNA-seq Kraken2 OR sylph OR "species-level" metagenomic profiling` |
| PubMed (E-utilities esearch, relevance sort) | curl | (1) `"RNA-seq"[Title/Abstract] AND "quality control"[Title/Abstract] AND (strandedness OR rRNA OR "gene body coverage" OR "transcript integrity")`; (2) `"sample swap"[Title/Abstract] AND "RNA-seq"[Title/Abstract]`; (3) `"sample identity"[Title/Abstract] AND "RNA-seq"[Title/Abstract] AND (genotype OR SNP)`; (4) `"RSeQC"[Title/Abstract] OR "transcript integrity number"[Title/Abstract]`; (5) `"Salmon"[Title/Abstract] AND ("run-to-run" OR nondeterministic OR non-deterministic OR "random seed") AND RNA-seq`; (6) `"MultiQC"[Title/Abstract] OR "nf-core"[Title/Abstract]`; (7) `"rRNA" AND "RNA-seq" AND "quality control" AND "threshold"` |
| Europe PMC (REST, published) | curl | 20 queries, e.g. `RSeQC comprehensive quality control RNA-seq experiments`; `transcript integrity number RNA-seq RSeQC TIN`; `Qualimap 2 RNA-seq gene body coverage quality`; `strandedness inference RNA-seq library type tool`; `ribosomal RNA contamination RNA-seq quality control depletion`; `RNA-seq sample swap identity genotype concordance detection`; `NGSCheckMate RNA-seq sample identity`; `sex check RNA-seq XIST Y chromosome sample mislabel`; `Kraken2 contamination RNA-seq species`; `sylph taxonomic profiling containment ANI`; `nf-core pipelines community curated reproducible`; `Salmon quasi-mapping non-deterministic reproducibility quantification`; `bioinformatics reproducibility containers software versions RNA-seq pipeline`; `MultiQC aggregate results bioinformatics analyses`; `fastp ultra-fast all-in-one FASTQ preprocessor`; `duplicate reads RNA-seq library complexity preseq`; `RNA-seq quality control guidelines outlier samples replicate correlation`; `SEQC MAQC-III RNA-seq quality assessment reproducibility`; `Sequencing quality control consortium RNA-seq standards ENCODE guidelines replicates`; `BLAS threads numerical reproducibility R bioinformatics`. Full list is in the scratch script, not in the repo. |
| bioRxiv | The bioRxiv MCP and API have no keyword search (date/category only), so preprints were searched through Europe PMC with `SRC:PPR` (which indexes bioRxiv/medRxiv/Research Square) | `(SRC:PPR) AND "RNA-seq" AND "quality control" AND (strandedness OR "gene body coverage" OR rRNA)`; `(SRC:PPR) AND RustQC`; `(SRC:PPR) AND "sample swap" AND "RNA-seq"`; `(SRC:PPR) AND Salmon AND reproducibility AND RNA-seq quantification`; `(SRC:PPR) AND "contamination" AND "RNA-seq" AND (Kraken2 OR sylph)`; `(SRC:PPR) AND "RNA-seq" AND "replicate" AND "correlation" AND "outlier"` |
| Tool documentation | direct fetch | RSeQC site (rseqc.sourceforge.net), Salmon docs (salmon.readthedocs.io), Picard metric definitions, dupRadar vignette, Kraken2 wiki manual, sylph docs, nf-core/rnaseq `usage.md`, `nextflow.config`, `nextflow_schema.json`, `CHANGELOG.md`, `docs/output.md`, ENCODE long-RNA data standards, RustQC README, MultiQC README, OpenBLAS USAGE/FAQ |
| GitHub | raw files and API | LICENSE files at tags (Salmon v1.0.0, v1.10.0, v1.10.3, v2.0.0), repo licence/latest-release for 25 repos, Salmon issues #185 and #562 |
| Package indexes (extra) | API | PyPI and bioconda for licence/version of 15+ tools |

Inclusion rule: a source is included if it states what a QC check measures, a tool, a threshold, a licence, or an evidence-based reproducibility practice. Exclusion: off-topic hits (non-RNA-seq, clinical case reports, single-cell-only, tissue-specific biology), duplicate records, and records with no QC content in the abstract.

## 2. Search log

| Stage | Count |
|---|---|
| Identified (titles returned and inspected): Europe PMC 112, PubMed E-utilities and Europe PMC preprints 51, PubMed MCP 21 | 184 |
| Additional tool/vendor pages and repos opened (documentation, LICENSE, releases, package indexes) | about 45 |
| Screened (title plus abstract or page text) | 184 records and about 45 pages |
| Read in depth (abstract or primary page quoted) | 33 |
| Included (cited below) | 33 |
| Excluded | 151 records |

Reported hit counts per query were much larger (for example 4,220 for the rRNA query), because I only inspected the top 6 to 20 per query. This is not an exhaustive systematic search; ranking was by relevance.

Exclusion reasons (counts approximate):
- Off-topic biology or method that merely mentions Qualimap/MultiQC/Kraken2 in methods (about 95): for example PMIDs 36555486, 35099752, 42380194, 42172220.
- Metagenome or microbiome application papers with no QC guidance for human RNA-seq (about 20).
- Clinical or ethics discussion of sex check (PMID 38822420, not a method; see section 6).
- DNA-methylation array identity QC (PMID 29881472, wrong data type).
- Proteogenomic mislabel challenge (PMID 34036290, wrong data type but relevant concept; excluded from the table).
- Conference abstracts and posters (about 8).
- Generic reproducibility papers with no RNA-seq-specific recommendation (about 10; two kept, 39123217 and 42030305).
- Preprints duplicating a later paper (PPR1045814 duplicates PMID 41593502).

Failed or unusable pages (listed so they can be retried):
- Europe PMC returned HTTP 503 on the first two runs. The cause was the `sort=RELEVANCE` parameter combined with a default user agent. Fixed by removing `sort` and sending a User-Agent. No data lost.
- GitHub API returned 403 "rate limit exceeded" partway through (RSEM repo, issue #562 body, Salmon issue comments via API). Salmon issue #185 was read through WebFetch instead. RSEM licence not checked.
- GitHub 404: `ewels/dupRadar`, `bioinformatics-centre/sylph`, `subread/featureCounts`, `rseqc/RSeQC`, `MolecularBiologyLaboratory/RSeQC`, `Teichlab/RSeQC`, `kundajilab/ataqv`, `Shicheng-Guo/ngscheckmate`. `ewels/MultiQC` and `ewels/RustQC` returned "Moved Permanently"; used `MultiQC/MultiQC` and `seqeralabs/RustQC`.
- Kraken2 wiki manual does not state a licence on that page. Sylph docs landing page: no contamination guidance, no RNA-seq statement, no licence statement.
- `pdftotext` is not installed, so the STAR manual PDF was downloaded but not read. STAR determinism is unverified.
- Salmon docs page: no mention of `lib_format_counts.json` or a deterministic/seed option found in the page text I searched (search of `seed`, `determin`, `reproduc` returned only the "randomized order" note).
- Qualimap site: version history and release date not retrieved. The Bitbucket LICENSE file was readable.
- The PubMed MCP translates natural language into noisy MeSH queries and returned unrelated items (for example PMIDs 36555486, 38822420).
- The PubMed query for Salmon non-determinism returned 0 papers: no peer-reviewed paper found on this topic in the search.
- No evidence found for RustQC in any preprint (0 hits); it is documented only on GitHub.

## 3. Findings per check

Licences below come from the tool's LICENSE file or package index on 2026-10-23. Where two sources disagree it is noted.

### 3.1 Strandedness inference

- **What it measures.** The orientation of reads relative to gene annotation, which sets the correct count and quantification setting.
- **Tools.** RSeQC `infer_experiment.py`; Salmon `-l A` auto-detection; `how_are_we_stranded_here`; Qualimap rnaseq; Picard `PCT_CORRECT_STRAND_READS`; resolveS (preprint).
- **Threshold stated.** nf-core/rnaseq defaults `stranded_threshold = 0.8` and `unstranded_threshold = 0.1` (V, `nextflow.config`). Quote (V, `usage.md`): "Forward stranded: At least 80% of the fragments are in the 'forward' orientation. Unstranded: The forward and reverse fractions differ by less than 10%. Undetermined: Samples that do not meet either criterion, possibly indicating issues such as genomic DNA contamination." Schema text (V): `stranded_threshold` is "must be at least 0.5". "These thresholds apply to both the strandedness inferred from Salmon outputs ... and how strandedness is inferred from RSeQC results."
- **nf-core cross-check (V, `docs/output.md`).** "If RSeQC disagrees with the expected strandedness, or returns 'undetermined' (which can indicate gDNA contamination), the summary row is marked as a fail." This is the practice our pipeline lacks: we rely on Salmon alone and have no second method.
- **RSeQC default** (V): `infer_experiment.py` runs on any BAM with a gene model BED; mapq default 30.
- **Source on cost of error** (V abstract, PMID 35065593): "when unknown or incorrectly specified can have detrimental effects on the reproducibility and accuracy of downstream analyses"; values outside the normal range "may indicate sample contamination". No numeric threshold in the abstract.
- **Licences.** RSeQC GPL-3.0-or-later (V, PyPI classifier, bioconda, site text "GNU General Public License (GPLv3)"); last release 5.0.6, 2026-10-19 on PyPI (V). Site release history says v5.0.5 on Jul. 26, 2026 (V). Salmon BSD-3-Clause from v2.0.0 (V). how_are_we_stranded_here MIT, v1.0.1 released 2021-03-15 (V). Picard MIT, 3.5.0 on 2026-07-31 (V). resolveS: preprint doi 10.64898/2026.01.08.698333, no licence or release checked (S).

### 3.2 rRNA content

- **Measures.** Fraction of reads or bases from ribosomal RNA; indicates failed depletion or poly-A selection.
- **Tools.** Picard `CollectRnaSeqMetrics` (`PCT_RIBOSOMAL_BASES`), RSeQC `read_distribution` plus the rRNA BED files it provides for human and mouse, SortMeRNA (nf-core default `ribo_removal_tool`, off by default: `remove_ribo_rna = false`, V), RNA-SeQC.
- **Threshold stated.** Picard definition (V): "Fraction of PF_ALIGNED_BASES that mapped to regions encoding ribosomal RNA"; threshold none. ENCODE long-RNA standards: no rRNA threshold (V, page stated none). The QC-DR paper (PMID 41593502, V abstract) found "% rRNA reads" among the most correlated metrics with sample quality but gives no cutoff in the abstract; full text not read (S for any number).
- **RSeQC page (V):** "We only provide rRNA bed files for human and mouse ... NO WARRANTY in completeness."
- **Licence.** Picard MIT (V); SortMeRNA licence not checked.

### 3.3 3' bias, TIN, RNA integrity

- **Measures.** Coverage skew along transcripts caused by degradation or oligo-dT priming.
- **Tools.** RSeQC `tin.py`; Picard `MEDIAN_3PRIME_BIAS`, `MEDIAN_5PRIME_TO_3PRIME_BIAS`, `MEDIAN_CV_COVERAGE`.
- **Threshold stated.** RSeQC TIN: "score (0 <= TIN <= 100) for each expressed transcript" and a sample-level median (medTIN) (V, RSeQC page). The page and the TIN abstract (PMID 26842848, V) give no pass/fail cutoff. Picard: "MEDIAN_CV_COVERAGE ... Ideal value = 0"; 3' bias is "mean coverage of the 3 prime-most 100 bases divided by the mean coverage of the whole transcript" computed over the 1000 most highly expressed transcripts; no thresholds (V).
- **Use of TIN in analysis.** The TIN paper (V abstract) reports that calibrating counts with TIN reduced false positives (sensitivity 0.89, specificity 0.91, accuracy 0.90 on spike-ins). This is an adjustment, not a gate.
- **Related evidence (V abstracts).** PMID 39715490: RIN correlated with 3' bias (r = -0.366 in autopsy cohorts). PMID 28446590: small RIN differences changed 8.1% of genes.
- **nf-core.** `tin` is available but not in the default `rseqc_modules` list (V schema: default `bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication`).
- **Licence.** RSeQC GPL-3+ (V); Picard MIT (V).

### 3.4 Gene-body coverage

- **Measures.** Uniformity of read coverage over scaled gene bodies.
- **Tools.** RSeQC `geneBody_coverage.py` (module listed in the RSeQC manual; S for its exact options), Qualimap rnaseq, RustQC (`rustqc rna`, V README lists Qualimap rnaseq reimplementation with "Gene body coverage").
- **Threshold stated.** None found. QC-DR (V abstract) introduced "Area Under the Gene Body Coverage Curve (AUC-GBC)" as a metric that correlated with quality, no cutoff stated in abstract.
- **Licence.** Qualimap GPL-2 (V, Bitbucket `LICENSE` text "GNU GENERAL PUBLIC LICENSE Version 2, June 1991"; bioconda says GPL-2.0-or-later, so uncertain whether "only" or "or later"); bioconda version 2.3 (V). Qualimap paper PMID 26428292 and 22914218 (V).

### 3.5 Duplication and library complexity

- **Measures.** PCR artefacts versus natural duplication of highly expressed genes; unique molecule count.
- **Tools.** dupRadar (R/Bioconductor, GPL-3, bioconda 1.40.0, V); Picard `MarkDuplicates` (`PERCENT_DUPLICATION`, `ESTIMATED_LIBRARY_SIZE`); RSeQC `read_duplication`; preseq (GPL-3.0, v3.2.0 released 2022-09-05, V, bioconda 3.2.0); RustQC reimplementations.
- **Threshold stated.** dupRadar vignette (V): helper lines "at the threshold of 1 read/bp and at 0.5 RPKM"; quote "The common denominator for both designs is the importance of having a low intercept"; the vignette "gives no numeric cutoff" for a good intercept or slope. Picard: none.
- **Pitfall (V, vignette).** In RNA-seq, high duplication at high expression is expected; a global duplicate fraction alone is misleading, which is why dupRadar stratifies by expression. dupRadar paper PMID 27769170 (V abstract).
- **Depth metric (V abstract, PMID 33712853).** MEND reads (mapped, exonic, non-duplicate) proposed for reporting depth; across 2,179 tumours: duplicates "3-100% of mapped reads (median 27%)", MEND "0-79% of total reads (median 50%)". These are observed distributions, not thresholds.
- **nf-core defaults (V).** `skip_preseq = true`, `skip_dupradar = false`, `skip_markduplicates = false`.
- **Licence.** dupRadar GPL-3 (V); preseq GPL-3 (V); Picard MIT (V).

### 3.6 Contamination (microbial, cross-sample)

- **Microbial.** Kraken2: default `--confidence` 0.0 (V quote "confidence score threshold (default: 0.0)"); standard database 29 GB as of Jan 2018 (V). Licence: repo LICENSE is MIT (V), bioconda metadata says GPL-3.0-or-later (conflict, S which applies to the packaged build); latest GitHub release v2.18.0-RC1 on 2026-10-20 (V; a release candidate), bioconda 2.17.2 (V). Sylph: repo MIT (V), bioconda 1.0.0 on 2026-09-14 (V). Sylph docs page states nothing about RNA-seq or thresholds (V). Benchmark (V abstract, PMID 42685266, 2026): for host-dominated RNA-seq, "assembly-based methods yield accurate taxonomic predictions but struggle at low read depth, whereas assembly-free methods are more robust in sparse settings at the cost of reduced precision", and a hybrid is proposed. No threshold. Implication: DNA-metagenome tools applied to human RNA-seq are research-grade, not standard QC.
- **Cross-sample.** NGSTroubleFinder (PMID 41608734, 2026): detects "cross-sample contamination in human whole-genome and whole-transcriptome sequencing data, sample swaps, and mismatches between the reported and the inferred genetic and transcriptomic sexes"; repo MIT, last push 2025-12-03 (V); runs on BAM/CRAM. IDCheck (PMID 23559639) also flags contamination from RNA-seq plus genotype.
- **Intrinsic signal.** RSeQC "undetermined" strandedness may flag gDNA (V nf-core text above).

### 3.7 Sex check

- **What sources say.** NGSTroubleFinder infers "genetic and transcriptomic sexes" from BAM and compares with reported sex (V abstract). No XIST or Y-gene thresholds are stated in the abstract. No source I opened gives an XIST/RPS4Y1 cutoff (threshold: none found).
- **Ethics note (V, PMID 38822420).** Sex-check QC "outs" transgender, gender diverse and intersex people and reinforces binary sex; the paper proposes alternative QC and clearer fields on requisition forms. This supports our decision to keep the sex check flag-only and to describe it as "genetic sex inferred from expression versus recorded sex" rather than a pass/fail.
- **Licence.** NGSTroubleFinder MIT (V).

### 3.8 Sample swap and identity

- **Tools.** NGSCheckMate (MIT, v1.0.1 on 2023-08-15, bioconda 1.0.1 on 2026-02, V). V README quotes: handles RNA-seq, FASTQ/BAM/VCF input, "alignment-free module ... less than one minute for RNA-seq using a single core", "Currently, it works only for human data". Paper (V abstract, PMID 28369524): "minimal requirement for sequencing depth (>0.5X)". IDCheck (PMID 23559639) needs genotype data. SMaSH (PMID 31888490) is also listed in PubMed (abstract not read). NGSTroubleFinder (MIT, above).
- **Requirement.** Any identity check needs an independent reference (genotype array, DNA-seq, or other samples from the same subject). Without it, only within-cohort duplicate/kinship detection is possible.
- **Threshold stated.** Not retrieved from the NGSCheckMate README beyond the depth figure (the correlation cutoff is in the paper's full text; not read).

### 3.9 Replicate correlation and cohort outliers

- **ENCODE (V, long-RNA data standards).** Quote: "Spearman correlation of >0.9 between isogenic replicates and >0.8 between anisogenic replicates". Also "Each replicate should have 30 million aligned reads" and "two or more replicates". No thresholds for mapping, rRNA, duplication, or strandedness.
- **Integration advice (V abstract, PMID 41593502).** "any individual QC metric is limited in its predictive value and suggests approaches based on the integration of multiple metrics with QC thresholds." Among the best correlated: % and number of uniquely aligned reads, % rRNA, number of detected genes, AUC-GBC. This supports our MAD-flag approach on several metrics at once, but numeric cutoffs are not in the abstract.
- **nf-core hard gate (V schema).** `min_mapped_reads` default 5: "Minimum percentage of uniquely mapped reads below which samples are removed from further processing."

### 3.10 Our existing metrics

- fastp: MIT, v1.4.0 released 2026-10-07 (V GitHub); fastp 1.0 paper PMID 41112039 (V title).
- No source I opened sets a Q30 or STAR unique % cutoff beyond nf-core's 5% mapping gate. Thresholds in our pipeline should be labelled project-defined.

### 3.11 Aggregation

- MultiQC: GPL-3-or-later (V README: "released under the GPL v3 or later licence"), v1.35 on 2026-05-13 (V). It can render RSeQC, Qualimap, Picard, Preseq, Salmon, featureCounts, dupRadar output (V nf-core `output.md`). MultiQC paper PMID 27312411; nf-core paper PMID 32055031 (V titles).
- RustQC: GPL-3.0+, v0.2.1 released 2026-04-09 (V). README claim: "All outputs are format- and numerically identical to the upstream tools" and "~14m 54s vs ~15h 34m" for a 186M-read BAM (V quotes; vendor benchmark, not independently verified). In nf-core/rnaseq, `--use_rustqc` is marked experimental (V schema). It is young (v0.2.x); treat as unproven.

## 4. Reproducibility practices

| Practice | What sources say | V/S |
|---|---|---|
| Containers | PMID 42030305 (2026): containers support reproducibility but "inappropriate or superficial use ... can undermine these benefits, leading to brittle environments, security risks, or false confidence in reproducibility." No detail of the nine tips read beyond the abstract. PMID 38025398 (CWL): combined with containerization, pipelines overcome "software incompatibility". | V (abstracts) |
| Version pinning, tests, validation | PMID 42756310 (2026): recommends "systematic validation checks (environment, inputs, and runtime)" and tests with real and synthetic data in CI. PMID 39123217: "genomic reproducibility" defined as "consistent results across technical replicates" and best practices recommended; specifics not read. | V (abstracts) |
| nf-core practice | `CHANGELOG.md` (V): "Improve reproducibility with Conda" (PR #1565); "set `LC_COLLATE=C` for reproducible gene-level output ordering" (PR #1883); sort inputs "keeping their cache keys stable across `-resume`" (PR #1885); "Exclude non-deterministic qualimap MultiQC files from snapshots" (PR #1770); "use file-names-only snapshots for all non-deterministic indices (STAR, Salmon, Kallisto, RSEM, HISAT2)" (PR #1687); "Deterministic mapping for STAR aligner" (issue #396, old release). Take-away: nf-core itself does not expect byte-identical Salmon/STAR indices and Qualimap outputs. | V |
| Salmon non-determinism | Issue #185 (V via WebFetch, open since 2018-01-17): maintainer "there is a small, inherent, non-determinism in the program, even if random seeds are all set identically ... due to the streaming nature of the first phase of the algorithm"; "deviations are expected to be very small"; (2018-05-29) "independent threads pull reads ... from a lock-free concurrent queue ... This is the major source of stochasticity." No seed option implemented per that page. Issue #562 (title only, V): different TPM between salmon 0.8.1 and 1.3.0. Salmon docs (V): reads should be in random order relative to targets; `-p/--threads` are used for quasi-mapping, quantification and bootstrapping. No peer-reviewed paper on this found. | V |
| Thread/BLAS determinism | OpenBLAS `USAGE.md` (V) documents `export OPENBLAS_NUM_THREADS=1` (line 37), but it does not state that doing so gives bit-identical results; I did not find a source that does. The PubMed hit "Bit-reproducible parallel phylogenetic tree inference" (PMID 41614360) shows the concept but is about phylogenetics (title only). R/DESeq2/limma sensitivity to BLAS was not verified. | Setting verified, claim of determinism S |
| STAR | `nf-core` issue #396 title "Deterministic mapping for STAR aligner" (V title only). STAR manual not read. | S |
| Hash manifests | No source consulted recommends hashing Salmon quant outputs. Given #185, a sha256 over `quant.sf` will vary between reruns; hashing inputs, references, versions, and config is robust. This is my inference. | S |

## 5. Licence table

| Tool | Licence | Latest release (date) | URL | V/S |
|---|---|---|---|---|
| RSeQC | GPL-3.0-or-later | 5.0.6 (PyPI upload 2026-10-19) | https://rseqc.sourceforge.net/ , PyPI RSeQC | V |
| Picard | MIT | 3.5.0 (2026-07-31) | https://github.com/broadinstitute/picard | V |
| Qualimap | GPL-2 (Bitbucket LICENSE; bioconda: GPL-2.0-or-later) | 2.3 (bioconda; release date not retrieved) | http://qualimap.conesalab.org/ | V licence, S date |
| dupRadar | GPL-3 | Bioconductor 1.40.0 (bioconda, 2026-02) | https://github.com/ssayols/dupRadar | V |
| preseq | GPL-3.0 | v3.2.0 (2022-09-05) | https://github.com/smithlabcode/preseq | V |
| Kraken2 | MIT in repo LICENSE; bioconda lists GPL-3.0-or-later | v2.18.0-RC1 (2026-10-20); stable on bioconda 2.17.2 | https://github.com/DerrickWood/kraken2 | V, conflict noted |
| sylph | MIT | latest tag "latest" (2026-09-14); bioconda 1.0.0 | https://github.com/bluenote-1577/sylph | V |
| RustQC | GPL-3.0+ | v0.2.1 (2026-04-09) | https://github.com/seqeralabs/RustQC | V |
| MultiQC | GPL-3.0-or-later | v1.35 (2026-05-13) | https://github.com/MultiQC/MultiQC | V |
| Salmon | v1.x and earlier: GPL-3 (LICENSE at v1.0.0, v1.10.0, v1.10.3). v2.0.0 and later: BSD-3-Clause | v2.8.0 (2026-09-25) | https://github.com/COMBINE-lab/salmon | V |
| fastp | MIT | v1.4.0 (2026-10-07) | https://github.com/OpenGene/fastp | V |
| nf-core/rnaseq | MIT | 3.27.0 (2026-09-23) | https://github.com/nf-core/rnaseq | V |
| STAR | MIT in current repo LICENSE; bioconda lists "GPLv3" for its build | 2.7.11b (2024-01-26) | https://github.com/alexdobin/STAR | V, conflict noted |
| featureCounts (Subread) | GPL-3.0-only (bioconda) | 2.1.1 (bioconda, 2025-05) | http://subread.sourceforge.net | S (not checked at source) |
| NGSCheckMate | MIT | v1.0.1 (2023-08-15) | https://github.com/parklab/NGSCheckMate | V |
| NGSTroubleFinder | MIT | no release metadata checked; last push 2025-12-03 | https://github.com/STALICLA-RnD/NGSTroubleFinder | V licence |
| how_are_we_stranded_here | MIT | v1.0.1 (2021-03-15) | https://github.com/signalbash/how_are_we_stranded_here | V |
| HTSinfer | Apache-2.0 (per paper abstract); PyPI classifier Apache; PyPI shows 0.2.0 dated 2020-06-04, which looks stale | S | https://github.com/zavolanlab/htsinfer (URL S) | S |

## 6. Recommendation

### 6.1 Checks to add, in order

1. **Second strandedness estimate from alignments, with a pass/fail against the Salmon call.** Highest value for lowest cost. RSeQC `infer_experiment.py` on the STAR BAM, with the nf-core thresholds 0.8 and 0.1 (V). Mark samples "undetermined" as possible gDNA. Today a single wrong Salmon library type changes counts silently.
2. **rRNA %, 3' bias, gene-body coverage, TIN.** Use Picard `CollectRnaSeqMetrics` (MIT) for `PCT_RIBOSOMAL_BASES`, `MEDIAN_3PRIME_BIAS`, `MEDIAN_CV_COVERAGE`, plus RSeQC `tin.py` and `geneBody_coverage.py` if you need curves. No source gives cutoffs, so flag by cohort MAD like the rest of our flags, and say so. The QC-DR paper supports combining several metrics rather than relying on one.
3. **Duplication and complexity.** dupRadar (intercept and slope; no numeric cutoff in the source) and MEND-read depth (PMID 33712853). Add `preseq` only for low-input protocols (nf-core skips it by default).
4. **Sample identity and swap.** NGSCheckMate (MIT) if you have any matched genotype, DNA, or multiple samples per subject; otherwise use it to catch accidental duplicates inside the cohort. NGSTroubleFinder (MIT) bundles identity, contamination, and expression-derived sex. Needs BAM and human data.
5. **Replicate correlation and hard gate.** Where replicates exist, report Spearman per pair and compare with the ENCODE 0.9 / 0.8 figures (V). Add a minimum-mapped gate like nf-core's 5% as an option; label it as inherited from nf-core, not validated for our data.
6. **Sex check.** Keep flag-only. Reword the output as expression-inferred sex versus recorded sex, with no XIST or Y-gene threshold cited because no source gave one. Do not remove samples on this basis alone (PMID 38822420).
7. **Microbial contamination.** Lowest priority. Kraken2 or sylph on unmapped reads is exploratory for human RNA-seq (benchmark PMID 42685266 reports a precision/sensitivity trade-off). Use only if a project question requires it.
8. **Reproducibility changes.** (a) Pin tool versions with a lock file and use container digests instead of floating conda specs (PMID 42030305 warns against false confidence). (b) Record `OPENBLAS_NUM_THREADS`/thread counts and Salmon `-p` in the manifest. (c) Hash inputs, references, versions, and parameters with sha256; for Salmon outputs, compare with a numeric tolerance instead of byte hashes (issue #185). (d) Stamp the Salmon version: v1.x is GPL-3 and v2.x is BSD-3 and quantification differs across versions (issue #562).

### 6.2 Exact tools and licence compatibility with our Salmon use

You wrote "GPL-3 Salmon use". That is true only for Salmon 1.x. Salmon 2.0.0 and later (current 2.8.0) is BSD-3-Clause (V). Check which version our conda env resolves, and update the licence notice in our docs.

Our pipeline invokes each tool as a separate process on files; it does not link their code. Under GPL-2 and GPL-3 this is generally treated as aggregation, so mixing licences across separately executed tools is normally not a problem. That is my reading, not legal advice. Redistribution of a container bundling these tools requires carrying each tool's licence text and source offer.

| Tool | Licence | Call it as a separate process? | Concern |
|---|---|---|---|
| Picard | MIT | yes | none |
| RSeQC | GPL-3+ | yes | none for separate execution; do not import its Python modules into our code unless our code is GPL-3-compatible |
| dupRadar | GPL-3 | yes (R script) | same as RSeQC |
| preseq | GPL-3 | yes | none |
| Qualimap | GPL-2 (maybe "or later") | yes | not compatible with GPL-3 code if linked; fine as a separate process |
| MultiQC | GPL-3+ | yes | none |
| NGSCheckMate, NGSTroubleFinder, how_are_we_stranded_here, sylph | MIT | yes | none |
| Kraken2 | MIT (repo) or GPL-3+ (bioconda metadata) | yes | check which applies if redistributing |
| RustQC | GPL-3+ | yes | experimental; do not make it the only implementation |
| Salmon 2.x | BSD-3 | yes | none; Salmon 1.x GPL-3 also fine |

Practical stack to add first: RSeQC (`infer_experiment`, `tin`, `geneBody_coverage`), Picard `CollectRnaSeqMetrics`, dupRadar, NGSCheckMate (and NGSTroubleFinder if sex/contamination is needed), all aggregated with MultiQC. All are free software with licences that allow separate-process use alongside Salmon.

### 6.3 What remains unverified

- Numeric thresholds for rRNA %, 3' bias, TIN, gene-body coverage, duplication, XIST/Y-gene sex calls, and NGSCheckMate correlation cutoffs: none stated by the sources I read (full text of QC-DR PMID 41593502, RNA-SeQC PMID 22539670, RSeQC TIN and NGSCheckMate papers not read; they may contain cutoffs or cohort-level observations).
- Whether any source states bit-level BLAS/OpenBLAS/thread behaviour for R or Python steps (DESeq2, PCA). OpenBLAS docs only show the thread setting.
- STAR run-to-run determinism (manual not read; nf-core issue #396 title only).
- featureCounts licence and version (taken from bioconda only).
- Qualimap exact licence variant ("only" versus "or later") and last release date.
- Kraken2 licence discrepancy (MIT in repo vs GPL-3+ in bioconda metadata) and STAR (MIT in repo vs GPLv3 on bioconda).
- HTSinfer licence and release status (Apache per abstract; PyPI metadata stale).
- Whether Salmon 2.x changed numerical output relative to 1.x for our data (only issue #562 about older versions seen).
- RustQC claims of numerical identity with upstream tools (vendor README only).
- Performance claims: none of the benchmark numbers were reproduced.
- The search was a relevance-ranked sample (top 6 to 20 per query), not a complete systematic review; no search of Embase, Scopus, or Web of Science.
