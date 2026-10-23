"""The pipeline as a chart: ordered stages with their Nextflow process names, and a status per stage computed from the task counts of a run.
Pure functions (no I/O) so they can be tested without a pipeline."""

PER_SAMPLE = {'VERIFY_FASTQ', 'FASTP', 'STAR_ALIGN', 'SAMTOOLS_SORT_INDEX', 'SAMTOOLS_STATS', 'FEATURECOUNTS', 'SALMON_QUANT', 'SAMPLE_QC_FLAGS'}

# lane: 'main' = one box in the main line; 'A' = STAR branch, 'B' = Salmon branch (both start after trimming and join at "Merge counts")
# what = one plain sentence for users; tool = the program that does it
STAGES = [
    dict(key='verify', label='Verify FASTQ', lane='main', procs=['VERIFY_FASTQ'], tool='gzip test',
         what='Checks that each FASTQ file is complete and readable before any heavy step starts.'),
    dict(key='fastp', label='Trim reads (fastp)', lane='main', procs=['FASTP'], tool='fastp',
         what='Removes adapters and low-quality bases and reads, and writes a quality report for every sample.'),
    dict(key='star', label='Align (STAR)', lane='A', procs=['STAR_ALIGN'], tool='STAR',
         what='Maps every read pair to the human genome (GRCh38) to find which gene it comes from.'),
    dict(key='sort', label='Sort + index (samtools)', lane='A', procs=['SAMTOOLS_SORT_INDEX', 'SAMTOOLS_STATS'], tool='samtools',
         what='Sorts and indexes the alignments so they can be counted, and records alignment statistics.'),
    dict(key='fc', label='Count genes (featureCounts)', lane='A', procs=['FEATURECOUNTS'], tool='featureCounts',
         what='Counts the read pairs that fall on each gene. These counts are the main data for the analysis.'),
    dict(key='salmon', label='Quantify (Salmon)', lane='B', procs=['SALMON_QUANT'], tool='Salmon',
         what='A second, independent way to measure expression (without alignment). It cross-checks the gene counts and tells the library type.'),
    dict(key='qc', label='Sample QC', lane='main', procs=['SAMPLE_QC_FLAGS', 'MERGE_SAMPLE_QC', 'CHECK_QC'], tool='QC rules',
         what='Marks each sample PASS, WARN or FAIL from its mapping rates and read quality. FAIL samples are left out of the analysis.'),
    dict(key='merge', label='Merge counts', lane='main', procs=['MERGE_FEATURECOUNTS', 'SALMON_TXIMPORT'], tool='merge + tximport',
         what='Combines the per-sample results into one gene-by-sample table (and imports the Salmon estimates).'),
    dict(key='cohort', label='Cohort QC', lane='main', procs=['COHORT_QC'], tool='R',
         what='Looks at all samples together: PCA, sample correlation, library sizes and outlier flags.'),
    dict(key='de', label='Differential expression', lane='main', procs=['DESEQ2_PRIMARY', 'DISCOVER_DESEQ2_PRIMARY', 'DISCOVER_DESEQ2_SECONDARY'], tool='DESeq2',
         what='Finds the genes whose expression differs between the two groups, adjusting for batch when one is given.'),
    dict(key='path', label='Pathways', lane='main', procs=['PATHWAY_ORA', 'PATHWAY_GSEA'], tool='clusterProfiler, fgsea',
         what='Asks which biological pathways are enriched among the differentially expressed genes.'),
    dict(key='cand', label='Candidate rules', lane='main', procs=['CANDIDATE_RULES'], tool='DESeq2 refits',
         what='Keeps only genes that are significant, well expressed, change in the same direction on both routes, and stay significant when each sample is left out in turn.'),
    dict(key='model', label='Model + cutoff', lane='main', procs=['FIT_MODEL', 'PACKAGE_MODEL'], tool='LASSO (glmnet)',
         what='Picks a small gene panel and a score cutoff from the candidates. This needs enough samples; with too few it is skipped with a clear message.'),
    dict(key='report', label='Report', lane='main', procs=['REPORT'], tool='report script',
         what='Collects the tables and figures into one HTML report with the QC findings and the limitations.'),
]


def compute(by_process, n_samples, job_status):
    """by_process: {PROCESS: {'done': n, 'active': n, 'failed': n}}. Returns one dict per stage with status pending|running|done|failed."""
    out = []
    for st in STAGES:
        done = sum(by_process.get(p, {}).get('done', 0) for p in st['procs'])
        active = sum(by_process.get(p, {}).get('active', 0) for p in st['procs'])
        failed = sum(by_process.get(p, {}).get('failed', 0) for p in st['procs'])
        exp = sum((n_samples if p in PER_SAMPLE else 1) for p in st['procs']) if n_samples else None
        if job_status == 'succeeded':
            status = 'done'
        elif failed:
            status = 'failed'
        elif exp is not None and done >= exp:
            status = 'done'
        elif active or done:
            status = 'running'
        else:
            status = 'pending'
        out.append(dict(key=st['key'], label=st['label'], lane=st['lane'], status=status, done=done, expected=exp, active=active, what=st['what'], tool=st['tool']))
    return out
