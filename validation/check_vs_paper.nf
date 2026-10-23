// NEW. Compares our results with the expected genes/pathways of Hosseini et al. 2023 AFTER the fact; never changes cutoffs.
process CHECK_VS_PAPER {
    label 'process_low'
    conda "${params.base_dir}/env/prefix/r_dge"
    publishDir "${params.outdir}/06_dge", mode: 'copy'

    input:
    path results_all
    path ora, stageAs: 'ora/*'
    path gsea, stageAs: 'gsea/*'

    output:
    path "check_vs_paper_genes.tsv", emit: genes
    path "check_vs_paper_pathways.tsv", emit: pathways
    path "check_vs_paper_summary.txt", emit: summary

    script:
    """
    check_vs_paper.R ${results_all} ora gsea ${params.de_padj}
    """

    stub:
    """
    printf "gene\\tmatch\\n" > check_vs_paper_genes.tsv; printf "pathway\\tmatch\\n" > check_vs_paper_pathways.tsv; echo stub > check_vs_paper_summary.txt
    """
}
