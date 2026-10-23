// NEW. Per-sample gate: waits until THIS sample's two FASTQs exist (download_run.sh only renames .part -> .fastq.gz
// after gzip -t passes), re-tests them, then releases the sample. Samples never wait for each other.
process VERIFY_FASTQ {
    tag "${meta.id}"
    label 'process_low'
    executor 'local'
    errorStrategy 'terminate'

    input:
    tuple val(meta), val(r1), val(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: reads

    script:
    """
    deadline=\$(( \$(date +%s) + ${params.wait_hours} * 3600 ))
    for f in ${r1} ${r2}; do
        until [ -s "\$f" ]; do
            if [ \$(date +%s) -ge \$deadline ]; then echo "ERROR: timed out waiting for \$f" >&2; exit 1; fi
            sleep 30
        done
        gzip -t "\$f"
    done
    ln -s ${r1} ${meta.id}_R1.fastq.gz
    ln -s ${r2} ${meta.id}_R2.fastq.gz
    """

    stub:
    """
    touch ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz
    """
}
