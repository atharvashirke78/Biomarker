// NEW (clusterProfiler / fgsea / msigdbr; written from package documentation, not copied).
process PATHWAY_ORA {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_pathway"
    publishDir "${params.outdir}/07_pathway/ora", mode: 'copy'

    input:
    path results_all

    output:
    path "ora_*.tsv", emit: tables
    path "*.png", optional: true, emit: plots
    path "ora_kegg_source.txt", emit: kegg_source

    script:
    """
    export R_USER_CACHE_DIR=${params.ref_dir}/r_cache; mkdir -p \$R_USER_CACHE_DIR   # msigdbr data cache lives with the reference, not on the small root disk
    ora.R ${results_all} ${params.de_padj} ${params.de_lfc} ${params.path_min_gs} ${params.path_max_gs} ${params.path_q} ${params.path_ora_plot_n} ${params.run_kegg}
    """

    stub:
    """
    for d in up down; do for s in GO_BP KEGG; do printf "ID\\tDescription\\tp.adjust\\n" > ora_\${d}_\${s}.tsv; done; done
    echo "stub" > ora_kegg_source.txt
    """
}

process PATHWAY_GSEA {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_pathway"
    publishDir "${params.outdir}/07_pathway/gsea", mode: 'copy'

    input:
    path results_all

    output:
    path "gsea_*.tsv", emit: tables
    path "*.png", optional: true, emit: plots

    script:
    """
    export R_USER_CACHE_DIR=${params.ref_dir}/r_cache; mkdir -p \$R_USER_CACHE_DIR
    gsea.R ${results_all} ${params.path_min_gs} ${params.path_max_gs} ${params.gsea_nperm} ${params.path_q} ${params.gsea_seed} ${params.path_gsea_plot_n} ${params.run_kegg}
    """

    stub:
    """
    for s in Hallmark KEGG GO_BP; do printf "pathway\\tNES\\tpadj\\n" > gsea_\${s}.tsv; done
    """
}
