// CHANGED: publish dir is a param (data/subsampled by default; test_data for the mini test); otherwise as written 07:06
// (first N read pairs via head; full files untouched).
process SUBSAMPLE {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.subsample_dir}", mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.sub.fastq.gz"), path("${meta.id}_R2.sub.fastq.gz"), emit: reads

    script:
    """
    set +o pipefail
    zcat ${r1} | head -n ${params.subsample_pairs * 4} | gzip > ${meta.id}_R1.sub.fastq.gz
    zcat ${r2} | head -n ${params.subsample_pairs * 4} | gzip > ${meta.id}_R2.sub.fastq.gz
    set -o pipefail
    n1=\$(zcat ${meta.id}_R1.sub.fastq.gz | wc -l); n2=\$(zcat ${meta.id}_R2.sub.fastq.gz | wc -l)
    if [ "\$n1" != "\$n2" ]; then echo "ERROR: subsample line counts differ \$n1 vs \$n2" >&2; exit 1; fi
    """

    stub:
    """
    touch ${meta.id}_R1.sub.fastq.gz ${meta.id}_R2.sub.fastq.gz
    """
}
