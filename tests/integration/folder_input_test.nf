// Test of modules/folder_input.nf only (no processes, no tools): scans --input_dir, prints the result, and checks that
// params.samplesheet can be set from it (QUANT_STAGE reads params.samplesheet).
// usage: nextflow run tests/integration/folder_input_test.nf -c nextflow/nextflow.config --code_dir <checkout> --input_dir <d> --output_dir <o> [--mode discover|apply] [-resume]
include { folderInput } from '../../modules/folder_input'

workflow {
    def r = folderInput(params.mode)
    params.samplesheet = r.sheet.toString()
    println "FOLDER_INPUT_OK sheet=${params.samplesheet} samples=${r.summary.n_samples} groups=${r.summary.groups} notes=${r.summary.notes}"
}
