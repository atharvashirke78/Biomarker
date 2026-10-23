// REPLACED (copy of align.nf STAR_ALIGN from /workspace/biomarker/pipeline had two-pass + sorted output; spec needs
// one-pass, unsorted BAM -> samtools sort, --outFilterMultimapNmax 20). STAR_INDEX takes sjdbOverhang from a fastp JSON.
process STAR_INDEX {
    label 'process_star_index'
    conda "${params.install_dir}/env/prefix/star"
    publishDir "${params.ref_dir}", mode: 'copy'

    input:
    path genome
    path gtf
    path fastp_json

    output:
    path "star_index", emit: index
    path "star_index_overhang.txt"

    script:
    """
    len=\$(grep -o '"read1_mean_length":[0-9]*' ${fastp_json} | head -1 | cut -d: -f2)
    ovh=\$(( len - 1 ))
    echo "read1_mean_length=\$len sjdbOverhang=\$ovh (from ${fastp_json})" | tee star_index_overhang.txt
    zcat ${gtf} > anno.gtf
    mkdir star_index
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star_index \\
        --genomeFastaFiles ${genome} --sjdbGTFfile anno.gtf --sjdbOverhang \$ovh \\
        --limitGenomeGenerateRAM ${params.star_index_ram}
    rm anno.gtf
    """

    stub:
    """
    mkdir star_index && touch star_index/Genome && echo "stub" > star_index_overhang.txt
    """
}

process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_star'
    conda "${params.install_dir}/env/prefix/star"
    publishDir "${params.outdir}/02_alignment", mode: 'copy', pattern: '*.Log.final.out',
        saveAs: { fn -> "${meta.run}/Log.final.out" }

    input:
    tuple val(meta), path(r1), path(r2)
    path index

    output:
    tuple val(meta), path("${meta.id}.Aligned.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}.Log.final.out"), emit: log

    script:
    """
    STAR --runThreadN ${task.cpus} --genomeDir ${index} \\
        --readFilesIn ${r1} ${r2} --readFilesCommand zcat \\
        --outFilterMultimapNmax ${params.star_multimap_max} \\
        --outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM \\
        --outFileNamePrefix ${meta.id}.
    """

    stub:
    """
    touch ${meta.id}.Aligned.out.bam
    printf "                          Uniquely mapped reads %% |\\t90.00%%\\n" > ${meta.id}.Log.final.out
    """
}
