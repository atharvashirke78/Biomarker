// NEW (written from samtools documentation, not copied).
process SAMTOOLS_SORT_INDEX {
    tag "${meta.id}"
    label 'process_samtools'
    conda "${params.install_dir}/env/prefix/samtools"
    publishDir "${params.outdir}/02_alignment", mode: 'copy', saveAs: { fn -> "${meta.run}/" + fn.replace(meta.id, meta.run) }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -m 2G -T tmp_${meta.id} -o ${meta.id}.sorted.bam ${bam}
    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    """
}

process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/samtools"
    publishDir "${params.outdir}/02_alignment", mode: 'copy', saveAs: { fn -> "${meta.run}/" + fn.replace(meta.id, meta.run) }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.flagstat"), path("${meta.id}.stats"), emit: stats

    script:
    """
    samtools flagstat ${bam} > ${meta.id}.flagstat
    samtools stats ${bam} > ${meta.id}.stats
    """

    stub:
    """
    touch ${meta.id}.flagstat ${meta.id}.stats
    """
}
