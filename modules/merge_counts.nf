// NEW. Cohort count matrices. Only samples that are not FAIL reach these processes (filtered in main.nf before collect()).
process MERGE_FEATURECOUNTS {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/04_counts", mode: 'copy'

    input:
    path counts, stageAs: 'fc/*'
    path sample_map

    output:
    path "counts_featurecounts.tsv", emit: matrix

    script:
    """
    merge_featurecounts.py fc counts_featurecounts.tsv ${sample_map}
    """

    stub:
    """
    # stub matrix uses the sample ids of the files this task was given (sorted), so a sample excluded by QC is missing here too and the contract/validation steps can be tested with -stub-run
    ls fc | sed 's/\\.featureCounts\\.txt\$//' | sort | awk 'BEGIN{printf "gene_id"} {n[NR]=\$1; printf "\\t%s", \$1} END{printf "\\n"; for (g=1; g<=2; g++) { printf "GENE%d", g; for (i=1;i<=NR;i++) printf "\\t%d", 40+g*10+i; printf "\\n" } }' > counts_featurecounts.tsv
    """
}

process SALMON_TXIMPORT {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_dge"
    publishDir "${params.outdir}/04_counts", mode: 'copy'

    input:
    path quant_dirs, stageAs: 'q/*'
    path tx2gene
    path sample_map

    output:
    path "counts_salmon_tximport.tsv", emit: counts
    path "salmon_tpm.tsv", emit: tpm
    path "salmon_tximport.rds", emit: rds

    script:
    // countsFromAbundance = "no" (spec); DESeq2 uses the rds via DESeqDataSetFromTximport
    """
    tximport.R ${tx2gene} q ${sample_map}
    """

    stub:
    """
    ls q | sed 's/_quant\$//' | sort | awk 'BEGIN{printf "gene_id"} {printf "\\t%s", \$1} END{printf "\\n"; for (g=1; g<=2; g++) { printf "GENE%d", g; for (i=1;i<=NR;i++) printf "\\t%.2f", 40.5+g*10+i; printf "\\n" } }' > counts_salmon_tximport.tsv
    cp counts_salmon_tximport.tsv salmon_tpm.tsv
    touch salmon_tximport.rds
    """
}
