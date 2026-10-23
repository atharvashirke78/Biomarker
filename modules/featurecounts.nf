// CHANGED vs the copy of /workspace/biomarker/pipeline/modules/featurecounts.nf: -Q 10, strandedness per sample from
// the same sample's Salmon lib_format_counts.json (ISF->-s 1, ISR->-s 2, else 0), gz GTF unpacked, outputs per sample to 04_counts,
// merge moved to merge_counts.nf.
process FEATURECOUNTS {
    tag "${meta.id}"
    label 'process_featurecounts'
    conda "${params.install_dir}/env/prefix/subread"
    publishDir "${params.outdir}/04_counts", mode: 'copy', saveAs: { fn -> "${meta.run}/" + fn.replace(meta.id, meta.run) }

    input:
    tuple val(meta), path(bam), path(bai), path(quant)
    path gtf

    output:
    tuple val(meta), path("${meta.id}.featureCounts.txt"), path("${meta.id}.featureCounts.txt.summary"), emit: counts

    script:
    """
    fmt=\$(grep -o '"expected_format": *"[A-Z]*"' ${quant}/lib_format_counts.json | sed -E 's/.*"([A-Z]+)"\$/\\1/')
    [ -n "\$fmt" ] || { echo "ERROR: no Salmon library format in ${quant}/lib_format_counts.json: strandedness cannot be inferred" >&2; exit 1; }
    case "\$fmt" in
        ISF|OSF|MSF) s=1 ;;
        ISR|OSR|MSR) s=2 ;;
        *)           s=0 ;;
    esac
    echo "salmon library format \$fmt -> featureCounts -s \$s" | tee ${meta.id}.strand.txt
    zcat -f ${gtf} > anno.gtf
    featureCounts -p --countReadPairs -B -C -Q ${params.fc_min_mapq} -s \$s -T ${task.cpus} \\
        -t exon -g gene_id -a anno.gtf -o ${meta.id}.featureCounts.txt ${bam}
    rm anno.gtf
    """

    stub:
    """
    printf "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${meta.id}.bam\\nGENE1\\t1\\t1\\t10\\t+\\t10\\t100\\nGENE2\\t1\\t1\\t10\\t+\\t10\\t50\\n" > ${meta.id}.featureCounts.txt
    printf "Status\\t${meta.id}.bam\\nAssigned\\t150\\nUnassigned_NoFeatures\\t10\\n" > ${meta.id}.featureCounts.txt.summary
    """
}
