// CHANGED: one process called once per results folder (+ once over everything); report goes to <folder>/qc/.
process MULTIQC {
    tag "${folder}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/multiqc"
    publishDir "${params.outdir}/${folder}/qc", mode: 'copy', pattern: 'multiqc_report.html'
    publishDir "${params.outdir}/${folder}", mode: 'copy', pattern: 'multiqc_report.html', saveAs: { fn -> folder == '08_multiqc' ? fn : null }

    input:
    tuple val(folder), path(qc_files, stageAs: 'qc/*')
    path sample_map

    output:
    tuple val(folder), path("multiqc_report.html"), emit: report

    script:
    """
    # sample names in the report = SRR (map: label and its tool-specific variants -> SRR)
    awk -F'\\t' -v OFS='\\t' 'NR>1 { print \$2, \$1; print \$2 "_quant", \$1; print \$2 ".sorted", \$1; print \$2 ".fastp", \$1; print \$2 ".featureCounts", \$1 }' ${sample_map} | awk -F'\\t' '{ print length(\$1) "\\t" \$0 }' | sort -k1,1nr | cut -f2- > replace_names.tsv   # longest names first (MultiQC replaces substrings in order)
    multiqc qc -n multiqc_report.html -f --title "${folder}" --replace-names replace_names.tsv
    """

    stub:
    """
    touch multiqc_report.html
    """
}

// QC rule: fastp html+json of ALL samples are copied into <each results folder>/qc/fastp/<SRR>.fastp.{html,json}
process QC_COPY_FASTP {
    label 'process_low'
    executor 'local'
    publishDir "${params.outdir}/01_fastp/qc", mode: 'copy'
    publishDir "${params.outdir}/02_alignment/qc", mode: 'copy'
    publishDir "${params.outdir}/03_quant_salmon/qc", mode: 'copy'
    publishDir "${params.outdir}/04_counts/qc", mode: 'copy'
    publishDir "${params.outdir}/05_cohort_qc/qc", mode: 'copy'
    publishDir "${params.outdir}/06_dge/qc", mode: 'copy'
    publishDir "${params.outdir}/07_pathway/qc", mode: 'copy'
    publishDir "${params.outdir}/08_multiqc/qc", mode: 'copy'

    input:
    path files, stageAs: 'in/*'
    path sample_map

    output:
    path "fastp", emit: dir

    script:
    """
    mkdir fastp
    tail -n +2 ${sample_map} | while IFS=\$'\\t' read -r srr label study group; do
        cp -L in/\${label}.fastp.html fastp/\${srr}.fastp.html
        cp -L in/\${label}.fastp.json fastp/\${srr}.fastp.json
    done
    """

    stub:
    """
    mkdir fastp
    tail -n +2 ${sample_map} | while IFS=\$'\\t' read -r srr label study group; do
        cp -L in/\${label}.fastp.html fastp/\${srr}.fastp.html
        cp -L in/\${label}.fastp.json fastp/\${srr}.fastp.json
    done
    """
}
