// Internal validation of the DISCOVER diagnostic model (upgrade H3c-H3f). ADDITIVE and REPORT-ONLY: reads the finished model, writes extra files,
// never feeds back into the model. Script: nextflow/bin/internal_validation.R (re-implements the recipe of fit_diagnostic_model.R inside every resample).
// LIMITATION (written into every output): the candidate screening (rules 1-3b) stays outside the resampling loop; NOT independent validation.
// Never fails the run: the script exits 0 and writes a 'NOT RUN: reason' row when there is no model, fewer than 10 samples or a class below 5;
// if R itself fails, the shell fallback writes the same three files, so the report (which waits for this step) is never blocked.
process INTERNAL_VALIDATION {
    label 'process_medium'
    conda "${params.install_dir}/env/prefix/r_classifier"
    publishDir "${params.outdir}${params.discover_dir}/04_validation", mode: 'copy'

    input:
    path counts
    path sheet
    path candidates
    path genes
    path spec
    path scores

    output:
    path "internal_validation_summary.tsv", emit: summary
    path "selection_frequency.tsv", emit: freq
    path "INTERNAL_VALIDATION.md", emit: md
    path "INTERNAL_VALIDATION.provenance.txt"

    script:
    """
    internal_validation.R --counts ${counts} --sheet ${sheet} --candidates ${candidates} --genes ${genes} --spec ${spec} --scores ${scores} \\
        --out . --repeats ${params.iv_repeats} --bootstrap ${params.iv_bootstrap} --subsamples ${params.iv_subsamples} || {
        r="NOT RUN: the internal-validation script failed (see .command.log); the model and all other outputs are unaffected"
        printf 'metric\\tvalue\\tci_low\\tci_high\\tlabel\\ninternal_validation\\tNA\\tNA\\tNA\\t%s\\n' "\$r" > internal_validation_summary.tsv
        printf 'gene_id\\tsymbol\\tdirection\\tfrequency\\tin_final_model\\n%s\\tNA\\tNA\\tNA\\tNA\\n' "\$r" > selection_frequency.tsv
        printf '# Internal validation of the diagnostic model\\n\\n**%s**\\n\\nThis is NOT independent validation.\\n' "\$r" > INTERNAL_VALIDATION.md; }
    { echo "## command (.command.sh)"; cat .command.sh
      echo "## sha256 inputs"; sha256sum ${counts} ${sheet} ${candidates} ${genes} ${spec} ${scores}
      echo "## sha256 outputs"; sha256sum internal_validation_summary.tsv selection_frequency.tsv INTERNAL_VALIDATION.md; } > INTERNAL_VALIDATION.provenance.txt
    """

    stub:
    """
    printf 'metric\\tvalue\\tci_low\\tci_high\\tlabel\\ninternal_validation\\tNA\\tNA\\tNA\\tNOT RUN: stub\\n' > internal_validation_summary.tsv
    printf 'gene_id\\tsymbol\\tdirection\\tfrequency\\tin_final_model\\n' > selection_frequency.tsv
    echo 'stub' > INTERNAL_VALIDATION.md
    touch INTERNAL_VALIDATION.provenance.txt
    """
}
