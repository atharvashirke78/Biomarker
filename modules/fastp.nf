// CHANGED vs 2026-10-08 07:06 version: spec cutoffs as params, 4 threads, per-sample publish layout,
// outputs fastp.html/json named <sample>.fastp.* (QC rule), meta-keyed json channel.
process FASTP {
    tag "${meta.id}"
    label 'process_fastp'
    conda "${params.install_dir}/env/prefix/fastp"
    publishDir "${params.outdir}/01_fastp", mode: 'copy', pattern: '*.fastp.{html,json}',
        saveAs: { fn -> "${meta.run}/" + fn.replace("${meta.id}.", '') }   // folder = SRR run accession (naming rule); label stays in meta.id

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.trimmed_1.fastq.gz"), path("${meta.id}.trimmed_2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json
    path "${meta.id}.fastp.json", emit: json_only
    path "${meta.id}.fastp.html", emit: html

    script:
    """
    fastp -i ${r1} -I ${r2} \\
        -o ${meta.id}.trimmed_1.fastq.gz -O ${meta.id}.trimmed_2.fastq.gz \\
        --detect_adapter_for_pe \\
        -q ${params.fastp_qual} --unqualified_percent_limit ${params.fastp_unqual_pct} \\
        --cut_tail --cut_window_size ${params.fastp_cut_window} --cut_mean_quality ${params.fastp_cut_mean_q} \\
        --length_required ${params.fastp_min_len} --n_base_limit ${params.fastp_n_base} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html 2> ${meta.id}.fastp.log
    """

    stub:
    """
    touch ${meta.id}.trimmed_1.fastq.gz ${meta.id}.trimmed_2.fastq.gz ${meta.id}.fastp.html
    echo '{"summary":{"before_filtering":{"read1_mean_length":100,"q30_rate":0.95,"total_reads":10},"after_filtering":{"q30_rate":0.96}},"filtering_result":{"passed_filter_reads":9}}' > ${meta.id}.fastp.json
    """
}
