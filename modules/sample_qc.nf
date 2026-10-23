// NEW. Per-sample QC flags (PASS/WARN/FAIL) from STAR, featureCounts, Salmon and fastp; stdlib python only.
process SAMPLE_QC_FLAGS {
    tag "${meta.id}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"

    input:
    tuple val(meta), path(star_log), path(fc_summary), path(quant), path(fastp_json)

    output:
    tuple val(meta), path("${meta.id}.qc.tsv"), emit: row

    script:
    """
    sample_qc.py '${meta.id}' '${meta.study}' '${meta.group}' ${star_log} ${fc_summary} ${quant} ${fastp_json} \\
        ${params.qc_star_warn} ${params.qc_fc_warn} ${params.qc_salmon_warn} ${params.qc_star_fail} ${params.qc_salmon_fail} > ${meta.id}.qc.tsv
    """

    stub:
    """
    printf "${meta.id}\\t${meta.study}\\t${meta.group}\\t90\\t80\\t80\\t95\\t90\\t${meta.id == params.stub_fail_sample ? 'FAIL' : 'PASS'}\\t${meta.id == params.stub_fail_sample ? 'stub test: marked FAIL' : ''}\\n" > ${meta.id}.qc.tsv
    """
}

process MERGE_SAMPLE_QC {
    label 'process_low'
    publishDir "${params.outdir}/04_counts", mode: 'copy'

    input:
    path rows, stageAs: 'rows/*'
    path sample_map

    output:
    path "sample_qc.tsv", emit: table

    script:
    // per-sample rows are keyed by label; the cohort table is keyed by SRR with label, study, group as extra columns
    """
    set -o pipefail
    export LC_ALL=C
    printf "srr\\tlabel\\tstudy\\tgroup\\tstar_unique_pct\\tfc_assigned_pct\\tsalmon_mapped_pct\\tfastp_q30_pct\\tfastp_pass_pct\\tstatus\\treason\\n" > sample_qc.tsv
    cat rows/*.qc.tsv | awk -F'\\t' -v OFS='\\t' 'NR==FNR{ if (FNR>1) srr[\$2]=\$1; next } { if (!(\$1 in srr)) { print "ERROR: no SRR for label " \$1 > "/dev/stderr"; exit 1 } print srr[\$1], \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10 }' ${sample_map} - | sort -k1,1 >> sample_qc.tsv
    """

    stub:
    """
    printf "srr\\tlabel\\tstudy\\tgroup\\tstar_unique_pct\\tfc_assigned_pct\\tsalmon_mapped_pct\\tfastp_q30_pct\\tfastp_pass_pct\\tstatus\\treason\\n" > sample_qc.tsv
    cat rows/*.qc.tsv | awk -F'\\t' -v OFS='\\t' 'NR==FNR{ if (FNR>1) srr[\$2]=\$1; next } { print srr[\$1], \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10 }' ${sample_map} - | sort -k1,1 >> sample_qc.tsv
    """
}
