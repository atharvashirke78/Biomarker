// OPT-IN (params.run_nested_validation, default false: it is slow). Fully nested internal validation, REPORT-ONLY: the candidate screen (rules 1-2) is
// repeated inside every training fold, plus a label-permutation null. Script: nextflow/bin/nested_full_validation.R. Never feeds back into the model.
// LIMITATION (in every output): rules 3a and 3b are not run, so it estimates the simplified procedure; NOT independent validation.
// Never fails the run: analysis problems give a 'NOT RUN: reason' row (exit 0); if R itself fails the shell fallback writes the same files.
process NESTED_VALIDATION {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_classifier"
    publishDir "${params.outdir}${params.discover_dir}/04_validation", mode: 'copy'

    input:
    path counts
    path sheet
    path rules
    path thresholds
    path questions

    output:
    path "nested_validation_summary.tsv", emit: summary
    path "nested_selection_frequency.tsv", emit: freq
    path "nested_permutation_null.tsv", emit: nulltab
    path "NESTED_VALIDATION.md", emit: md

    script:
    """
    nested_full_validation.R --counts ${counts} --sheet ${sheet} --questions ${questions} --rules ${rules} --thresholds ${thresholds} \\
        --repeats ${params.nested_repeats} --permutations ${params.nested_permutations} --outer-folds ${params.nested_outer_folds} --cores ${task.cpus < 4 ? task.cpus : 4} \\
        --rscript-dge ${params.install_dir}/env/prefix/r_dge/bin/Rscript --rscript-clf ${params.install_dir}/env/prefix/r_classifier/bin/Rscript --out . || {
        r="NOT RUN: the nested-validation script failed (see .command.log); the model and all other outputs are unaffected"
        printf 'metric\\tvalue\\tci_low\\tci_high\\tlabel\\nnested_validation\\tNA\\tNA\\tNA\\t%s\\n' "\$r" > nested_validation_summary.tsv
        printf 'gene_id\\tfrequency\\tlabel\\n' > nested_selection_frequency.tsv
        printf 'permutation\\tauc\\tlabel\\n' > nested_permutation_null.tsv
        printf '# Nested internal validation\\n\\n**%s**\\n\\nThis is NOT independent validation.\\n' "\$r" > NESTED_VALIDATION.md; }
    """

    stub:
    """
    printf 'metric\\tvalue\\tci_low\\tci_high\\tlabel\\nnested_validation\\tNA\\tNA\\tNA\\tNOT RUN: stub\\n' > nested_validation_summary.tsv
    printf 'gene_id\\tfrequency\\tlabel\\n' > nested_selection_frequency.tsv
    printf 'permutation\\tauc\\tlabel\\n' > nested_permutation_null.tsv
    echo stub > NESTED_VALIDATION.md
    """
}
