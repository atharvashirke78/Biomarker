// -entry RUN --mode apply: score the QUANT counts with the LOCKED diagnostic model (lock/diagnostic_v1). Nothing is fitted or tuned.
// Wraps scripts/apply_locked_model.R unchanged. The script reads and hash-checks "lock/diagnostic_v1" relative to its working folder, so the
// lock folder of this checkout is staged into the task folder under that exact name. QC FAIL samples are excluded by the script itself.
// Outputs: <outdir>/09_locked_model/scores/{scores.tsv, excluded_samples.tsv, run_info.tsv [, evaluation.tsv]}.
process LOCKED_SCORE {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/r_classifier"
    publishDir "${params.outdir}/09_locked_model", mode: 'copy'

    input:
    path counts
    path sheet, stageAs: 'sample_sheet_apply.tsv'
    path qc, stageAs: 'qc_table.tsv'
    path lock, stageAs: 'lock/diagnostic_v1'
    val cohort            // development | validation (written into every output; development is labelled optimistic)
    val evaluate          // true only when every sample has a known label and both groups are present

    output:
    path "scores", emit: scores

    script:
    def ev = evaluate ? '--evaluate' : ''
    """
    Rscript ${params.code_dir}/scripts/apply_locked_model.R --counts ${counts} --sample-sheet ${sheet} --qc-table ${qc} \\
        --cohort ${cohort} ${ev} --out-dir scores
    """

    stub:
    """
    mkdir scores
    printf 'sample_id\\tscore\\tcall\\n' > scores/scores.tsv
    """
}
