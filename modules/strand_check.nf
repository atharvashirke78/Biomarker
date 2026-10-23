// Independent strandedness check (flag only). featureCounts -s comes from the Salmon library type (modules/featurecounts.nf); here RSeQC
// infer_experiment.py reads the sorted BAM against a BED12 built from the same GTF and the two calls are compared. It never changes a count, a
// QC status or an exclusion; a disagreement is listed in 05_cohort_qc/strandedness_check.tsv and as a warning in the report.
process GTF_TO_BED12 {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"

    input:
    path gtf

    output:
    path "gene_model.bed12", emit: bed

    script:
    """
    gtf_to_bed12.py ${gtf} gene_model.bed12
    """

    stub:
    """
    printf "chr1\\t0\\t100\\tTX1\\t0\\t+\\t0\\t100\\t0\\t1\\t100,\\t0,\\n" > gene_model.bed12
    """
}

process STRAND_CHECK {
    tag "${meta.id}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/rseqc"

    input:
    tuple val(meta), path(bam), path(bai), path(quant)
    path bed

    output:
    tuple val(meta), path("${meta.id}.strand_check.tsv"), emit: row

    script:
    // python3 of the rseqc environment runs strand_check.py (standard library only); no other environment is needed
    """
    infer_experiment.py -r ${bed} -i ${bam} -s ${params.strand_sample_reads} > ${meta.id}.infer_experiment.txt
    strand_check.py '${meta.id}' ${meta.id}.infer_experiment.txt ${quant} ${params.strand_stranded_min} ${params.strand_unstranded_max_diff} > ${meta.id}.strand_check.tsv
    """

    stub:
    """
    printf "${meta.id}\\t0.9\\t0.05\\t0.05\\tforward\\tISF\\tforward\\t1\\tyes\\t\\n" > ${meta.id}.strand_check.tsv
    """
}

process MERGE_STRAND_CHECK {
    label 'process_low'
    publishDir "${params.outdir}/05_cohort_qc", mode: 'copy'

    input:
    path rows, stageAs: 'rows/*'
    path sample_map

    output:
    path "strandedness_check.tsv", emit: table

    script:
    // sample ids are labels in the per-sample rows; the cohort tables are keyed by SRR (same rule as MERGE_SAMPLE_QC)
    """
    set -o pipefail
    export LC_ALL=C
    printf "srr\\tlabel\\trseqc_forward\\trseqc_reverse\\trseqc_failed\\trseqc_call\\tsalmon_format\\tsalmon_call\\tfeaturecounts_s\\tagree\\tnote\\n" > strandedness_check.tsv
    cat rows/*.strand_check.tsv | awk -F'\\t' -v OFS='\\t' 'NR==FNR{ if (FNR>1) srr[\$2]=\$1; next } { s = (\$1 in srr) ? srr[\$1] : \$1; \$1 = s OFS \$1; print }' ${sample_map} - | sort -k1,1 >> strandedness_check.tsv
    """

    stub:
    """
    printf "srr\\tlabel\\trseqc_forward\\trseqc_reverse\\trseqc_failed\\trseqc_call\\tsalmon_format\\tsalmon_call\\tfeaturecounts_s\\tagree\\tnote\\n" > strandedness_check.tsv
    cat rows/*.strand_check.tsv >> strandedness_check.tsv
    """
}
