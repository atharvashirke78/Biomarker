// CHANGED: kept, but OFF by default (params.run_fastqc = false) because fastp already gives raw+trimmed QC.
process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/fastqc"
    publishDir "${params.outdir}/01_fastp/${meta.run}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${r1} ${r2}
    """

    stub:
    """
    for f in ${r1} ${r2}; do b=\$(basename \$f .fastq.gz); touch \${b}_fastqc.zip \${b}_fastqc.html; done
    """
}
