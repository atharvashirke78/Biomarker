// NEW (refactor Step 2; design docs/two_pipeline_design.md sections 2, 3, 2.5). Contract files, ANALYSIS input validation, manifests.
// Written for this project; stdlib python scripts in nextflow/bin.

// QUANT -> ANALYSIS contract: <outdir>/contract/{counts_matrix.tsv, counts_salmon_tximport.tsv, salmon_tpm.tsv, salmon_tximport.rds, sample_sheet.tsv, qc_table.tsv}.
// The 01_-08_ folders are unchanged; these are copies/derivations.
process MAKE_CONTRACT {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/contract", mode: 'copy'

    input:
    path resolved_sheet
    path sample_qc
    path cohort_flags
    path primary_matrix, stageAs: 'primary/*'      // featureCounts matrix (STAR route) or tximport counts (Salmon-only)
    path tximport_counts
    path tximport_tpm
    path tximport_rds
    val count_source

    output:
    path "counts_matrix.tsv", emit: counts
    path "counts_salmon_tximport.tsv", emit: sal_counts
    path "salmon_tpm.tsv"
    path "salmon_tximport.rds", emit: rds
    path "sample_sheet.tsv", emit: sheet
    path "qc_table.tsv", emit: qc

    script:
    """
    cp -L ${primary_matrix} counts_matrix.tsv
    [ "${tximport_counts}" = counts_salmon_tximport.tsv ] || cp -L ${tximport_counts} counts_salmon_tximport.tsv
    [ "${tximport_tpm}" = salmon_tpm.tsv ] || cp -L ${tximport_tpm} salmon_tpm.tsv
    [ "${tximport_rds}" = salmon_tximport.rds ] || cp -L ${tximport_rds} salmon_tximport.rds
    make_contract.py ${resolved_sheet} ${sample_qc} ${cohort_flags} counts_matrix.tsv ${count_source} sample_sheet.tsv qc_table.tsv
    """

    stub:
    """
    cp -L ${primary_matrix} counts_matrix.tsv
    [ "${tximport_counts}" = counts_salmon_tximport.tsv ] || cp -L ${tximport_counts} counts_salmon_tximport.tsv
    [ "${tximport_tpm}" = salmon_tpm.tsv ] || cp -L ${tximport_tpm} salmon_tpm.tsv
    [ "${tximport_rds}" = salmon_tximport.rds ] || cp -L ${tximport_rds} salmon_tximport.rds
    make_contract.py ${resolved_sheet} ${sample_qc} ${cohort_flags} counts_matrix.tsv ${count_source} sample_sheet.tsv qc_table.tsv
    """
}

// ANALYSIS input checks (design section 3). A failed check stops the run with the messages in .command.err / the Nextflow error;
// nothing downstream starts. Outputs the validated sample sheet that DESEQ2 reads.
process VALIDATE_ANALYSIS_INPUT {
    tag "${route}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/00_analysis_input", mode: 'copy', saveAs: { fn -> "${prefix}_${fn}" }

    input:
    val prefix
    path counts
    path sheet, stageAs: 'in_sample_sheet.tsv'      // distinct from the output name: never write through a staged symlink
    path symbols
    path qc_table
    val route
    val secondary

    output:
    path "sample_sheet.validated.tsv", emit: sheet
    path "analysis_input_check.json", emit: check

    script:
    def qc = qc_table.name == 'NO_FILE' ? '' : "--qc_table ${qc_table}"
    """
    validate_analysis_input.py --counts ${counts} --sheet ${sheet} --symbols ${symbols} --route ${route} --secondary ${secondary} \\
        --min_samples ${params.analysis_min_samples_per_group} --max_unmapped_pct ${params.analysis_max_unmapped_gene_pct} ${qc} \\
        --out_sheet sample_sheet.validated.tsv --out_json analysis_input_check.json
    """

    stub:
    """
    validate_analysis_input.py --counts ${counts} --sheet ${sheet} --symbols ${symbols} --route ${route} --secondary ${secondary} \\
        --min_samples ${params.analysis_min_samples_per_group} --max_unmapped_pct ${params.analysis_max_unmapped_gene_pct} \\
        --out_sheet sample_sheet.validated.tsv --out_json analysis_input_check.json
    """
}

// Manifest writer. run_meta / params_json are JSON strings made in main.nf (resolved parameters = every params.* value of this run).
process WRITE_MANIFEST {
    tag "${stage}"
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/contract", mode: 'copy'

    input:
    val stage
    val run_meta
    val params_json
    path files, stageAs: 'files/*'
    val file_roles           // list of "role=filename" matching the staged files
    val fastq                // list of absolute FASTQ paths (QUANT) or []
    path quant_manifest
    val with_refs

    output:
    path "manifest_${stage}.json", emit: manifest

    script:
    def roles = file_roles.collect { "--file ${it.split('=')[0]}=files/${it.split('=')[1]}" }.join(' ')
    def fq = fastq ? "--fastq " + fastq.join(' ') : ''
    def refs = with_refs ? "--refs ${params.ref_dir}" : ''
    def qm = quant_manifest.name == 'NO_FILE' ? '' : "--quant_manifest ${quant_manifest}"
    """
    cat > run_meta.json <<'META_EOF'
${run_meta}
META_EOF
    cat > params.json <<'PARAMS_EOF'
${params_json}
PARAMS_EOF
    write_manifest.py ${stage} run_meta.json params.json manifest_${stage}.json ${roles} ${fq} ${refs} ${qm} --base ${params.base_dir} --envbase ${params.install_dir}
    """

    stub:
    """
    cat > run_meta.json <<'META_EOF'
${run_meta}
META_EOF
    cat > params.json <<'PARAMS_EOF'
${params_json}
PARAMS_EOF
    python3 -c "import json; json.load(open('run_meta.json')); json.load(open('params.json'))"   # stub still checks that the resolved parameters are valid JSON
    echo '{"stub": true, "stage": "${stage}"}' > manifest_${stage}.json
    """
}
