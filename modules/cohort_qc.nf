// REPLACED qc_pca.nf (old version: study-based flags, Salmon only). New: PCA coloured by group / shaped by study,
// correlation heatmap, library-size barplot, 3-MAD flags from the group centre (flag only).
process COHORT_QC {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_dge"
    publishDir "${params.outdir}/05_cohort_qc", mode: 'copy'

    input:
    path counts
    path sample_map
    path sample_qc

    output:
    path "pca_group_study.png", emit: pca
    path "sample_correlation_heatmap.png", emit: heatmap
    path "library_size_barplot.png", emit: libsize
    path "cohort_qc_flags.tsv", emit: flags
    path "cohort_qc_summary.json", emit: summary
    path "cohort_qc_sex.tsv", emit: sex

    script:
    """
    cohort_qc.R ${counts} ${sample_map} ${sample_qc} ${params.cohort_min_count} ${params.cohort_min_samples} ${params.cohort_outlier_mad}
    """

    stub:
    """
    touch pca_group_study.png sample_correlation_heatmap.png library_size_barplot.png cohort_qc_flags.tsv cohort_qc_sex.tsv
    echo '{"stub": true}' > cohort_qc_summary.json
    """
}
