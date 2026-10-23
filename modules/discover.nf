// -entry DISCOVER (Phase 1, 2026-10-15): input contract -> DESeq2 (one or two routes) -> candidate rules 1-3b -> diagnostic model fit -> packaging.
// Scripts: nextflow/bin/{check_discover_input.py, deseq2.R (unchanged), candidate_rules.R, fit_diagnostic_model.R, package_model.py}.
// Every process writes <PROCESS>.provenance.txt: the exact command script (.command.sh), tool versions and sha256 of its inputs and outputs.
// Outputs go to params.outdir (a new run folder); nothing is written into lock/.

process DISCOVER_INPUT_CHECK {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}${params.discover_dir}/00_input", mode: 'copy'

    input:
    path counts
    path clinical
    path symbols
    val counts_secondary      // absolute path or '' (read, not staged: only its header is needed)
    val holdout               // absolute path or ''

    output:
    path "input_checks.json", emit: checks
    path "sample_sheet.discover.tsv", emit: sheet
    path "counts_training.tsv", emit: counts_training
    path "training_samples.txt"
    path "holdout_samples.txt"
    path "DISCOVER_INPUT_CHECK.provenance.txt"

    script:
    def sec = counts_secondary ? "--counts-secondary ${counts_secondary}" : ''
    def ho = holdout ? "--holdout ${holdout}" : ''
    """
    check_discover_input.py --question ${params.question} --counts ${counts} ${sec} --clinical ${clinical} ${ho} --symbols ${symbols} \\
        --min-samples ${params.analysis_min_samples_per_group} --max-unmapped-pct ${params.analysis_max_unmapped_gene_pct} \\
        --robustness \$(grep -E '^robustness_method:' ${params.questions_config} | awk '{print \$2}') \\
        --loo-max-samples \$(grep -E '^robustness_loo_max_samples:' ${params.questions_config} | awk '{print \$2}')
    { echo "## command (.command.sh)"; cat .command.sh; echo "## tool versions"; python3 --version; python3 -c 'import pandas; print("pandas", pandas.__version__)'
      echo "## sha256 inputs"; sha256sum ${counts} ${clinical} ${symbols} ${counts_secondary ?: ''} ${holdout ?: ''} ${params.questions_config}
      echo "## sha256 outputs"; sha256sum input_checks.json sample_sheet.discover.tsv counts_training.tsv training_samples.txt holdout_samples.txt; } > DISCOVER_INPUT_CHECK.provenance.txt
    """

    stub:
    """
    echo '{"stub": true}' > input_checks.json
    cp ${clinical} sample_sheet.discover.tsv
    cp ${counts} counts_training.tsv
    touch training_samples.txt holdout_samples.txt DISCOVER_INPUT_CHECK.provenance.txt
    """
}

process DISCOVER_DESEQ2 {
    tag "${prefix}"
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_dge"
    publishDir "${params.outdir}${params.discover_dir}/01_de", mode: 'copy'

    input:
    val prefix
    val mode
    path counts_tsv
    path counts_rds
    path sample_sheet
    path symbols

    output:
    path "${prefix}_results_all.tsv", emit: results
    path "${prefix}_summary.json", emit: summary
    path "${prefix}_*.tsv"
    path "${prefix}_*.png"
    path "${prefix}.provenance.txt"

    script:
    """
    COUNTS=${counts_tsv}; [ "${mode}" = tximport ] && COUNTS=${counts_rds}
    deseq2.R ${mode} \$COUNTS ${sample_sheet} ${symbols} ${prefix} \\
        ${params.de_min_count} ${params.de_min_samples} ${params.de_padj} ${params.de_lfc} ${params.de_volcano_labels} ${params.de_heatmap_genes} ${params.de_volcano_seed}
    { echo "## command (.command.sh)"; cat .command.sh; echo "## tool versions"
      Rscript -e 'cat(R.version.string, "| DESeq2", as.character(packageVersion("DESeq2")), "| apeglm", as.character(packageVersion("apeglm")), if (requireNamespace("tximport", quietly = TRUE)) paste("| tximport", packageVersion("tximport")), "\\n")'
      echo "## sha256 inputs"; sha256sum \$COUNTS ${sample_sheet} ${symbols}
      echo "## sha256 outputs"; sha256sum ${prefix}_results_all.tsv ${prefix}_summary.json ${prefix}_significant_up.tsv ${prefix}_significant_down.tsv ${prefix}_normalized_counts.tsv; } > ${prefix}.provenance.txt
    """

    stub:
    """
    touch ${prefix}_results_all.tsv ${prefix}_stub.png ${prefix}.provenance.txt
    echo '{"stub": true}' > ${prefix}_summary.json
    """
}

process CANDIDATE_RULES {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_dge"
    publishDir "${params.outdir}${params.discover_dir}/02_candidates", mode: 'copy'

    input:
    path counts
    path sheet
    path de_primary
    path de_secondary         // assets/NO_FILE when there is one route
    path rules
    path thresholds
    path questions

    output:
    path "candidate_table.tsv", emit: table
    path "candidate_rules_method_note.txt", emit: note
    path "funnel_counts.tsv"
    path "loo_pass_matrix.tsv"
    path "robustness_refits.tsv", emit: refits
    path "robustness_feasibility.tsv"
    path "reproduction_check_19_sample_fit.tsv"
    path "loo_failures_by_dropped_sample.tsv", optional: true
    path "sessionInfo.txt"
    path "tool_versions.txt"
    path "CANDIDATE_RULES.provenance.txt"

    script:
    def sec = de_secondary.name == 'NO_FILE' ? '' : "--de-secondary ${de_secondary}"
    """
    candidate_rules.R --counts ${counts} --sheet ${sheet} --de-primary ${de_primary} ${sec} --rules ${rules} --thresholds ${thresholds} --questions ${questions} --patient-level-cv ${params.patient_level_cv} --de-padj ${params.de_padj} --de-min-count ${params.de_min_count} --de-min-samples ${params.de_min_samples}
    { echo "## command (.command.sh)"; cat .command.sh; echo "## tool versions"; cat tool_versions.txt
      echo "## sha256 inputs"; sha256sum ${counts} ${sheet} ${de_primary} ${de_secondary.name == 'NO_FILE' ? '' : de_secondary} ${rules} ${thresholds} ${questions}
      echo "## sha256 outputs"; sha256sum candidate_table.tsv funnel_counts.tsv loo_pass_matrix.tsv reproduction_check_19_sample_fit.tsv candidate_rules_method_note.txt \$(ls loo_failures_by_dropped_sample.tsv 2>/dev/null); } > CANDIDATE_RULES.provenance.txt
    """

    stub:
    """
    touch candidate_table.tsv candidate_rules_method_note.txt funnel_counts.tsv loo_pass_matrix.tsv robustness_feasibility.tsv reproduction_check_19_sample_fit.tsv sessionInfo.txt tool_versions.txt CANDIDATE_RULES.provenance.txt robustness_refits.tsv
    """
}

process FIT_MODEL {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_classifier"
    publishDir "${params.outdir}${params.discover_dir}/03_model", mode: 'copy'

    input:
    path counts
    path sheet
    path candidates
    path de_primary
    path questions
    path method_note

    output:
    path "model_genes.tsv", emit: genes
    path "model_spec.txt", emit: spec
    path "inner_fold_cutoffs.tsv", emit: folds
    path "development_scores.tsv", emit: scores
    path "development_performance.tsv", emit: perf
    path "development_sens_spec.tsv", emit: sens_spec
    path "tool_versions.txt"
    path "FIT_MODEL.provenance.txt"

    script:
    """
    fit_diagnostic_model.R --counts ${counts} --sheet ${sheet} --candidates ${candidates} --de ${de_primary} --questions ${questions} --method-note ${method_note} --patient-level-cv ${params.patient_level_cv}
    { echo "## command (.command.sh)"; cat .command.sh; echo "## tool versions"; cat tool_versions.txt
      echo "## sha256 inputs"; sha256sum ${counts} ${sheet} ${candidates} ${de_primary} ${questions} ${method_note}
      echo "## sha256 outputs"; sha256sum model_genes.tsv model_spec.txt inner_fold_cutoffs.tsv development_scores.tsv development_performance.tsv; } > FIT_MODEL.provenance.txt
    """

    stub:
    """
    touch development_sens_spec.tsv model_genes.tsv model_spec.txt inner_fold_cutoffs.tsv development_scores.tsv development_performance.tsv tool_versions.txt FIT_MODEL.provenance.txt
    """
}

process PACKAGE_MODEL {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}${params.discover_dir}/03_model", mode: 'copy'

    input:
    path genes
    path spec
    path folds
    path scores
    path perf
    path checks
    path note
    path thresholds
    path questions
    path counts
    path clinical
    path candidates

    output:
    path "MODEL_MANIFEST.tsv"
    path "DISCOVER_SUMMARY.md"
    path "PACKAGE_MODEL.provenance.txt"

    script:
    """
    package_model.py --model-files ${genes} ${spec} ${folds} ${scores} ${perf} --input-checks ${checks} --method-note ${note} \\
        --thresholds ${thresholds} --questions ${questions} --run-name ${workflow.runName} \\
        --inputs counts=${counts} clinical=${clinical} candidate_table=${candidates} thresholds=${thresholds} questions=${questions}
    { echo "## command (.command.sh)"; cat .command.sh; echo "## tool versions"; python3 --version
      echo "## sha256 outputs"; sha256sum MODEL_MANIFEST.tsv DISCOVER_SUMMARY.md; } > PACKAGE_MODEL.provenance.txt
    """

    stub:
    """
    touch MODEL_MANIFEST.tsv DISCOVER_SUMMARY.md PACKAGE_MODEL.provenance.txt
    """
}
