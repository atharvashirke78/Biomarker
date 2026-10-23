// -entry RUN: final report. Wraps scripts/make_report.py unchanged. The script reads a run-folder layout (input/, pipeline/<stage>/..., run_info.json),
// so this process rebuilds that layout in its own task folder from STAGED inputs (never from the published output folder, which may still be
// being written) and runs the script there. Missing optional parts (e.g. discover files in apply mode) arrive as NO_FILE placeholders and are
// removed before the script runs.
// Outputs: <outdir>/10_report/{report/report.md, report/report.html, tables/, plots/, model/, manifest.json, tool_versions.tsv}; tool_versions.tsv also in <outdir>.
process REPORT {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/10_report", mode: 'copy'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'tool_versions.tsv'

    input:
    path run_info, stageAs: 'run_info.json'
    path input_files, stageAs: 'input/*'
    path qc_table, stageAs: 'pipeline/quant/04_counts/sample_qc.tsv'
    path de_files, stageAs: 'pipeline/discover/01_de/*'
    path cand_files, stageAs: 'pipeline/discover/02_candidates/*'
    path model_files, stageAs: 'pipeline/discover/03_model/*'
    path scores, stageAs: 'pipeline/apply/scores'
    path cohort_plots, stageAs: 'pipeline/quant/05_cohort_qc/*'
    path dge_plots, stageAs: 'pipeline/analysis/06_dge/*'
    path pathway_plots, stageAs: 'pipeline/analysis/07_pathway/*'
    path quant_manifest, stageAs: 'pipeline/quant/manifest_quant.json'
    path validation_files, stageAs: 'pipeline/discover/04_validation/*'      // internal validation files (upgrade H3), or the NO_FILE placeholder

    output:
    path "report", emit: report
    path "tables"
    path "plots"
    path "model", optional: true
    path "manifest.json"
    path "tool_versions.tsv"

    script:
    """
    find pipeline input -name 'NO_FILE*' -delete
    find pipeline -type d -empty -delete
    python3 -I ${params.code_dir}/scripts/make_report.py --run .
    """

    stub:
    """
    mkdir report tables plots
    echo '# stub report' > report/report.md
    touch manifest.json tool_versions.tsv
    """
}
