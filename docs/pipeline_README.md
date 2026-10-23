# Biomarker RNA-seq pipeline, diagnostic branch (FASTQ -> pathway analysis)

Nextflow DSL2. Entry: `nextflow/main.nf`, config: `nextflow/nextflow.config` (every cutoff is a documented param), modules: `modules/*.nf`, R/python helpers: `nextflow/bin/`, pinned conda envs: `env/*.yml`.
Run with the wrapper `env/bin/nextflow` (Java 17 + Nextflow 24.10.5 in a micromamba env; put `/workspace/Biomarker_pipeline/env/bin` on PATH).

## Input
`samplesheets/diagnostic_20.csv` (sample,run,study,group): 20 paired-end SRA runs, 10 tumour / 10 normal, studies SRP219837 and SRP301216.
FASTQs: `data/tumor_sample/` and `data/normal_sample/` as `<SRR>_R1.fastq.gz`, `<SRR>_R2.fastq.gz` (never modified).

## Reference (ref/, built once; skipped if present)
- genome: `ref/GRCh38_primary_chr.fa`, chr1-22,X,Y,M cut read-only from /workspace/genome/hg38/genome/Homo_sapiens_assembly38.fasta (headers reduced to contig names). Differs from GENCODE's primary-assembly genome only by hard-masked N (13.9 Mb, mostly chrY PARs, centromere-adjacent regions of chr5/14/19/21/22) and 94 IUPAC-vs-N positions; no real base differences (ref/check/).
- annotation: GENCODE v50 GTF restricted to chr1-22,X,Y,M, the 63 chrY genes lying entirely in PAR1/PAR2 removed (chrX copies carry PAR counts); transcript FASTA restricted to the transcript IDs of that GTF; tx2gene and gene_id->symbol from the same GTF. Counts of everything removed: docs/run_log.md.
- Salmon index: decoy-aware (decoys = the 25 genome contigs, gentrome = transcripts + genome, k=31, --gencode). STAR index: sjdbOverhang = read length - 1 from the fastp JSON of the first sample.

## Flow
Per sample (independent, no waiting for other samples): VERIFY_FASTQ -> [SUBSAMPLE] -> FASTP -> STAR_ALIGN -> SAMTOOLS_SORT_INDEX -> SAMTOOLS_STATS; FASTP -> SALMON_QUANT; BAM + Salmon lib type -> FEATURECOUNTS (strandedness per sample); STAR/featureCounts/Salmon/fastp -> SAMPLE_QC_FLAGS (PASS/WARN/FAIL; FAIL samples are excluded from the cohort steps).
Cohort (after collect()): MERGE_FEATURECOUNTS and SALMON_TXIMPORT -> COHORT_QC (PCA, correlation, library sizes, 3-MAD flags) -> DESEQ2 (design ~ study + group, reference normal; featureCounts route primary, Salmon route second + concordance) -> PATHWAY_ORA, PATHWAY_GSEA -> CHECK_VS_PAPER -> MULTIQC (8 reports + combined).

## Profiles
`-profile conda` (envs from env/*.yml), `-profile conda,stubtest` (`-stub-run`, own work/out dirs), `-profile conda,test` (4 samples x 100,000 pairs, work_test/results_test, relaxed DESeq2 pre-filter via params.test_mode), default = full 20-sample run into results/.

## Naming rule: SRR everywhere
All per-sample outputs and every cohort sample ID use the SRA run accession (SRRxxxxxxx = `meta.run`), never the label (T1_A, N1_A ...). `results/sample_map.tsv` (srr, label, study, group; generated from `samplesheets/diagnostic_20.csv` by MAKE_SAMPLE_MAP) links them. Per-sample folders/files: `01_fastp/<SRR>/fastp.html|json`, `02_alignment/<SRR>/<SRR>.sorted.bam`, `03_quant_salmon/<SRR>/`, `04_counts/<SRR>/<SRR>.featureCounts.txt`; QC copies `<folder>/qc/fastp/<SRR>.fastp.{html,json}`; cohort tables (counts matrices, tximport, DESeq2 colData and normalized counts, sample_qc.tsv with extra columns label/study/group, MultiQC sample names) are keyed by SRR. Implementation: per-sample processes keep the label as their internal key and rename only in `publishDir` (so the Nextflow cache/`-resume` is unaffected); cohort processes receive `sample_map.tsv`.
Resuming: `nextflow run ... -resume <session-id>` (see `nextflow log`); a bare `-resume` resumes the last session, which may be a test run.

## Results layout
`results/01_fastp ... 08_multiqc`, each with `qc/` (fastp html/json of all samples + that folder's MultiQC report), `results/RESULTS_INDEX.md`.
