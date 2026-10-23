// CHANGED: QUANT_SALMON renamed SALMON_QUANT (spec); publishes whole quant dir per sample to 03_quant_salmon/<sample>/;
// tximport moved to merge_counts.nf. Method unchanged from the copy of /workspace/biomarker/pipeline/modules/salmon.nf.
process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_salmon'
    conda "${params.install_dir}/env/prefix/salmon"
    publishDir "${params.outdir}/03_quant_salmon", mode: 'copy', saveAs: { fn -> fn.replaceFirst("^${meta.id}_quant", "${meta.run}") }

    input:
    tuple val(meta), path(r1), path(r2)
    path index

    output:
    tuple val(meta), path("${meta.id}_quant"), emit: quant

    script:
    """
    salmon quant -i ${index} -l A -1 ${r1} -2 ${r2} -p ${task.cpus} \\
        --validateMappings --gcBias --seqBias -o ${meta.id}_quant
    """

    stub:
    """
    mkdir -p ${meta.id}_quant/aux_info
    printf "Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\nENST1\\t1000\\t900\\t1\\t10\\n" > ${meta.id}_quant/quant.sf
    echo '{"num_processed": 10, "num_mapped": 8, "percent_mapped": 80.0}' > ${meta.id}_quant/aux_info/meta_info.json
    echo '{"expected_format": "IU", "compatible_fragment_ratio": 1.0}' > ${meta.id}_quant/lib_format_counts.json
    """
}
