// Extra RNA-QC metrics (flag only, cohort-relative; never a gate): RSeQC read distribution, duplication, TIN and gene-body coverage, plus rRNA and
// mitochondrial share from the featureCounts file. Slow tools run on a seeded random subset of transcripts (params.rnaqc_subset_transcripts,
// params.rnaqc_genebody_transcripts). Nothing here changes a count, a QC status or an exclusion; flags go to 05_cohort_qc/rna_qc_metrics.tsv and the report.
process RNAQC_SUBSET {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"

    input:
    path bed

    output:
    path "subset_tin.bed12", emit: tin
    path "subset_genebody.bed12", emit: genebody

    script:
    """
    rnaqc_subset_bed.py ${bed} subset_tin.bed12 ${params.rnaqc_subset_transcripts} ${params.rnaqc_seed}
    rnaqc_subset_bed.py ${bed} subset_genebody.bed12 ${params.rnaqc_genebody_transcripts} ${params.rnaqc_seed}
    """

    stub:
    """
    cp ${bed} subset_tin.bed12; cp ${bed} subset_genebody.bed12
    """
}

process RNAQC_SAMPLE {
    tag "${meta.id}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/rseqc"

    input:
    tuple val(meta), path(bam), path(bai), path(fc_counts)
    path bed
    path bed_tin
    path bed_genebody
    path symbols

    output:
    tuple val(meta), path("${meta.id}.rnaqc.tsv"), emit: row

    script:
    def gb = params.rnaqc_genebody ? "geneBody_coverage.py -i ${bam} -r ${bed_genebody} -o ${meta.id}.gb > ${meta.id}.gb.log 2>&1; GB=${meta.id}.gb.geneBodyCoverage.txt" : "GB=NONE"
    """
    export MPLBACKEND=Agg
    read_distribution.py -i ${bam} -r ${bed} > ${meta.id}.rd.txt
    read_duplication.py -i ${bam} -o ${meta.id}.dup > ${meta.id}.dup.log 2>&1
    tin.py -i ${bam} -r ${bed_tin} > ${meta.id}.tin.log 2>&1
    ${gb}
    rnaqc_row.py '${meta.id}' ${meta.id}.rd.txt ${meta.id}.dup.pos.DupRate.xls ${bam.baseName}.summary.txt \$GB ${fc_counts} ${symbols} > ${meta.id}.rnaqc.tsv
    """

    stub:
    """
    printf "${meta.id}\\t80.00\\t5.00\\t10.00\\t30.00\\t1.000\\t0.500\\t1.000\\t\\n" > ${meta.id}.rnaqc.tsv
    """
}

process MERGE_RNAQC {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/05_cohort_qc", mode: 'copy'

    input:
    path rows, stageAs: 'rows/*'
    path sample_map

    output:
    path "rna_qc_metrics.tsv", emit: table

    script:
    """
    cat rows/*.rnaqc.tsv > all_rows.tsv
    rnaqc_flags.py all_rows.tsv ${sample_map} ${params.cohort_outlier_mad} ${params.rnaqc_min_samples} > rna_qc_metrics.tsv
    """

    stub:
    """
    cat rows/*.rnaqc.tsv > all_rows.tsv
    rnaqc_flags.py all_rows.tsv ${sample_map} ${params.cohort_outlier_mad} ${params.rnaqc_min_samples} > rna_qc_metrics.tsv
    """
}
