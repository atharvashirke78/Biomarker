// Checks parseSamplesheet() + fastqPair() of nextflow/main.nf on a sheet without running any process (the QUANT memory guard is not
// involved: no workflow of main.nf runs). usage (project root):
//   env/bin/nextflow run tests/integration/parse_sheet_test.nf -c nextflow/nextflow.config --samplesheet <sheet.csv> [--legacy_samplesheet false]
include { parseSamplesheet; fastqPair } from '../../nextflow/main.nf'

workflow {
    def rows = parseSamplesheet(params.samplesheet)
    rows.each { r -> def p = fastqPair(r); println "ROW\t${r.sample}\t${r.group}\t${r.baseline}\t${r.batch}\t${p[0]}\t${p[1]}" }
    println "N_ROWS\t${rows.size()}"
}
