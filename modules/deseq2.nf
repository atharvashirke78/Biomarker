// REPLACED (old: tximport object only, no QC cutoffs). One process, called twice from main.nf (featureCounts route and
// Salmon route). design ~ batch + group, reference = baseline group of the sample sheet (refactor Step 2; was ~ study + group, reference normal).
// Written from DESeq2 vignette/documentation.
process DESEQ2 {
    tag "${prefix}"
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_dge"
    publishDir "${params.outdir}/06_dge", mode: 'copy'

    input:
    val prefix
    val mode
    path counts_tsv
    path counts_rds
    path sample_sheet
    path symbols

    output:
    path "${prefix}_results_all.tsv", emit: results
    path "${prefix}_significant_up.tsv", emit: up
    path "${prefix}_significant_down.tsv", emit: down
    path "${prefix}_normalized_counts.tsv"
    path "${prefix}_volcano.png", emit: volcano
    path "${prefix}_top50_heatmap.png", emit: heatmap
    path "${prefix}_MA.png", emit: ma
    path "${prefix}_summary.json", emit: summary

    script:
    """
    COUNTS=${counts_tsv}; [ "${mode}" = tximport ] && COUNTS=${counts_rds}
    deseq2.R ${mode} \$COUNTS ${sample_sheet} ${symbols} ${prefix} \\
        ${params.de_min_count} ${params.de_min_samples} ${params.de_padj} ${params.de_lfc} ${params.de_volcano_labels} ${params.de_heatmap_genes} ${params.de_volcano_seed}
    """

    stub:
    """
    for f in results_all significant_up significant_down normalized_counts; do printf "gene_id\\tgene_name\\n" > ${prefix}_\$f.tsv; done
    touch ${prefix}_volcano.png ${prefix}_top50_heatmap.png ${prefix}_MA.png
    echo '{"stub": true}' > ${prefix}_summary.json
    """
}

process DE_CONCORDANCE {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/r_dge"
    publishDir "${params.outdir}/06_dge", mode: 'copy'

    input:
    path primary
    path secondary

    output:
    path "concordance_featurecounts_vs_salmon.tsv", emit: table

    script:
    """
    concordance.R ${primary} ${secondary} ${params.de_padj} ${params.de_lfc}
    """

    stub:
    """
    printf "metric\\tvalue\\nstub\\t1\\n" > concordance_featurecounts_vs_salmon.tsv
    """
}
