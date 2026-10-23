// One-workflow entry (-entry RUN), step 1: input_dir (paired FASTQ + metadata.csv) -> samplesheet and clinical table, built automatically.
// Implemented as a Groovy function, not a process: QUANT_STAGE reads the samplesheet with Groovy at workflow-definition time
// (parseSamplesheet), so the sheet must exist before the stages are defined. The scan itself is scripts/input_scan.py (unchanged logic;
// opt-in flags --batch-optional and --min-per-group). Any input problem stops the run here with the scanner's message.
//
// Output (<output_dir>/00_input/): scan_files.tsv, problem_files.tsv, metadata.used.csv, metadata.resolved.tsv, group_mapping.tsv,
// samplesheet_quant.csv, clinical_discover.tsv, sample_sheet_apply.tsv, input_summary.json.
// Resume: the scan is repeated into a temporary folder and must equal the stored 00_input; if the folder changed since the run
// started, the run stops (start a new output_dir) so that a resumed run never mixes two inputs.

def folderInput(String mode) {
    if (!params.input_dir) error "-entry RUN needs --input_dir <folder with paired FASTQ files and metadata.csv>"
    if (!(mode in ['discover', 'apply'])) error "--mode ${mode}: -entry RUN supports discover or apply"
    def inDir = file(params.input_dir)
    if (!inDir.isDirectory()) error "--input_dir ${params.input_dir} is not a folder"
    if (!params.output_dir) error "no output folder: pass --output_dir (default is <input_dir>/Result)"
    def outDir = file(params.output_dir)
    outDir.mkdirs()
    def dest = file("${outDir}/00_input")
    def py = "${params.install_dir}/env/prefix/pyqc/bin/python3"
    def scanner = "${params.code_dir}/scripts/input_scan.py"
    [py, scanner].each { if (!file(it).exists()) error "input scan: ${it} not found (install_dir = ${params.install_dir}, code_dir = ${params.code_dir})" }
    if (dest.exists() && !workflow.resume) error "${dest} exists: this output_dir already holds a run. Use -resume to continue it, or choose a new --output_dir"

    def tmp = file("${outDir}/.input_scan_${System.currentTimeMillis()}")
    // discover scans as discover (patient_id needed, groups checked); apply on FASTQ scans as quant (metadata with a label per sample and the baseline group are needed to run QUANT)
    def scanMode = mode == 'discover' ? 'discover' : 'quant'
    def cmd = [py, '-I', scanner, '--input-dir', inDir.toString(), '--mode', scanMode, '--out-dir', tmp.toString(),
               '--config', "${params.code_dir}/conf/input.config", '--batch-optional']
    if (mode == 'discover') cmd += ['--patient-optional', '--min-per-group', "${params.analysis_min_samples_per_group}".toString(), '--questions', params.questions_config.toString()]
    def pr = new ProcessBuilder(cmd.collect { it.toString() }).redirectErrorStream(true).start()
    def out = pr.inputStream.text
    def rc = pr.waitFor()
    if (rc != 0) {
        error "INPUT CHECK FAILED (nothing was run). ${out.trim()}\n(scanner files, including any metadata.DRAFT.csv: ${tmp})"
    }

    def same = { f -> dest.resolve(f).exists() && dest.resolve(f).text == tmp.resolve(f).text }
    if (dest.exists()) {
        if (!['samplesheet_quant.csv', 'metadata.resolved.tsv', 'clinical_discover.tsv'].every(same))
            error "the input folder changed since this output_dir was started (${dest} differs from a new scan; scan kept in ${tmp}). Start a new --output_dir"
        tmp.deleteDir()
    } else {
        tmp.moveTo(dest)
    }
    log.info out.trim()
    def summary = new groovy.json.JsonSlurper().parse(file("${dest}/input_summary.json"))
    return [dir: dest, sheet: file("${dest}/samplesheet_quant.csv"), clinical: file("${dest}/clinical_discover.tsv"),
            apply_sheet: file("${dest}/sample_sheet_apply.tsv"), resolved: file("${dest}/metadata.resolved.tsv"), summary: summary]
}


// ---- run start: lock verification, code version, pod. Stops the run if the locked-model files do not match their hashes.
// Writes <output_dir>/run_info.json (one attempt appended per start / resume) and <output_dir>/params.json (all parameters of this launch).
def runStartInfo(String mode, String paramsText) {
    def sh = { List c ->
        def p = new ProcessBuilder(c.collect { it.toString() }).directory(new File(params.code_dir.toString())).redirectErrorStream(true).start()
        [p.inputStream.text.trim(), p.waitFor()]
    }
    def lock = sh(['sha256sum', '-c', 'lock/diagnostic_v1/sha256sums_final.txt'])
    def nOk = lock[0].readLines().count { it.endsWith(': OK') }
    def nAll = lock[0].readLines().findAll { it.contains(': ') }.size()
    if (lock[1] != 0) error "LOCK CHECK FAILED: the locked-model files do not match lock/diagnostic_v1/sha256sums_final.txt (${nOk} of ${nAll} OK). Nothing was run.\n${lock[0]}"
    def rd = { String p -> new File(p).exists() ? new File(p).text.trim() : null }
    def attempt = [kind: workflow.resume ? 'resume' : 'new', start: workflow.start?.toString(), stop_after: params.stop_after, status: 'started']
    def jf = new File("${params.output_dir}/run_info.json")
    def info
    if (jf.exists() && workflow.resume) {
        info = new groovy.json.JsonSlurper().parse(jf)
        info.attempts << attempt
    } else {
        def md = new File("${params.input_dir}/${params.metadata_file}")
        info = [run_id: workflow.runName, mode: mode, input_dir: params.input_dir.toString(), output_dir: params.output_dir.toString(), counts: null,
                command: workflow.commandLine, started: workflow.start?.toString(),
                git_commit: sh(['git', 'rev-parse', 'HEAD'])[0], git_branch: sh(['git', 'branch', '--show-current'])[0],
                uncommitted_files: sh(['git', 'status', '--porcelain'])[0].readLines(),
                lock_diagnostic_v1: "${nOk}/${nAll} OK".toString(), metadata_sha256: md.exists() ? sh(['sha256sum', md.toString()])[0].split(' ')[0] : null,
                pod: [nproc: Runtime.runtime.availableProcessors(), cgroup_memory_max: rd('/sys/fs/cgroup/memory.max'), cpu_max: rd('/sys/fs/cgroup/cpu.max')],
                label: 'development / technical output; not clinically validated', attempts: [attempt]]
    }
    jf.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
    new File("${params.output_dir}/params.json").text = paramsText
    return info
}


// clinical rows -> only samples that are columns of the QUANT counts matrix (QC FAIL samples are not in it)
process FILTER_CLINICAL {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/00_input", mode: 'copy'

    input:
    path clinical
    path counts

    output:
    path "clinical_discover.qcpass.tsv", emit: clinical
    path "excluded_after_quant.tsv", emit: excluded

    script:
    """
    filter_clinical.py --clinical ${clinical} --counts ${counts}
    """

    stub:
    """
    cp ${clinical} clinical_discover.qcpass.tsv
    head -1 ${clinical} > excluded_after_quant.tsv
    """
}


// after quantification: the samples excluded by the sample QC are listed (00_input/qc_excluded.tsv); the run stops with a clear message
// if fewer than min_total samples (cohort QC) or fewer than min_per_group per group (discover) are left
process CHECK_QC {
    label 'process_low'
    conda "${params.install_dir}/env/prefix/pyqc"
    publishDir "${params.outdir}/00_input", mode: 'copy'

    input:
    path sample_qc
    val min_per_group

    output:
    path "qc_excluded.tsv", emit: excluded

    script:
    """
    check_qc_min.py --qc ${sample_qc} --min-total 3 --min-per-group ${min_per_group}
    """

    stub:
    """
    check_qc_min.py --qc ${sample_qc} --min-total 3 --min-per-group ${min_per_group}
    """
}


// ---- startup checks (-entry RUN): everything the run needs, checked BEFORE any compute, printed in plain language, ALL problems at once.
// [ OK ] fine | [WARN] the run can start, but look at it | [MISSING] the run cannot work: it stops after the list.
// Under -stub-run the software / reference / memory checks are skipped (tests use fake references); --check_only true runs them anyway and stops.
def startupChecks(String mode) {
    def L = []
    def ok   = { String item, String d -> L << ['[ OK ]', item, d] }
    def warn = { String item, String d -> L << ['[WARN]', item, d] }
    def bad  = { String item, String d -> L << ['[MISSING]', item, d] }
    def gb   = { long b -> String.format('%.1f GB', b / 1e9) }

    // ---- input folder: the first thing that is checked, and the only thing if it is not set
    def inputOk = false
    if (!params.input_dir) {
        bad('input folder', 'NO INPUT FOLDER SET. Set input_dir in the USER SETTINGS block at the top of nextflow/nextflow.config, or start with --input_dir <folder>. The pipeline never uses a default folder or test data.')
    } else if (!new File(params.input_dir.toString()).isDirectory()) {
        bad('input folder', "${params.input_dir} does not exist or is not a folder (check the spelling and that the disk is mounted)")
    } else {
        def d = new File(params.input_dir.toString())
        def fq = d.listFiles().findAll { it.name ==~ /.*\.(fastq|fq)\.gz$/ }
        if (!fq) bad('FASTQ files', "no *.fastq.gz or *.fq.gz files in ${d}")
        else ok('FASTQ files', "${fq.size()} file(s) in ${d} (pairs are checked in the next step)")
        if (!new File(d, params.metadata_file.toString()).isFile()) warn(params.metadata_file.toString(), "not found in ${d}: the next step stops and writes a draft (metadata.DRAFT.csv) to fill in; every run needs the columns sample_id and group (see docs/RUN_IN_5_STEPS.md)")
        else ok(params.metadata_file.toString(), 'found')
        inputOk = true
        // output folder: writable (nearest existing parent), nothing is created here
        def o = new File((params.output_dir ?: "${d}/Result").toString()); def anc = o
        while (anc != null && !anc.exists()) anc = anc.parentFile
        if (anc == null || !anc.canWrite()) bad('output folder', "${o}: no write permission (nearest existing folder: ${anc})")
        else ok('output folder', "${o} (results, logs and the work folder go here)")
        if (anc != null && fq) {
            def nPairs = Math.max(1, (int) (fq.size() / 2)); def need = nPairs * 13_000_000_000L + 5_000_000_000L
            if (anc.usableSpace < need) warn('free disk space', "${gb(anc.usableSpace)} reported for ${anc}; about ${gb(need)} are needed (about 13 GB of work space per sample). This pod may report the shared volume, not your quota: check your quota. This is a warning, not a stop.")
            else ok('free disk space', "${gb(anc.usableSpace)} reported (about ${gb(need)} needed; the shared volume may report cluster space, not your quota)")
        }
    }

    if (params.check_only || !workflow.stubRun) {
        // ---- software
        def vOk = workflow.nextflow.version.matches('>=24.04.0')
        if (vOk) ok('nextflow', "version ${workflow.nextflow.version}") else bad('nextflow', "version ${workflow.nextflow.version} is too old: 24.04 or newer is needed")
        def inst = params.install_dir.toString()
        if (!new File("${inst}/env/bin/micromamba").exists()) bad('micromamba', "${inst}/env/bin/micromamba not found: the install folder is wrong or incomplete (install_dir = ${inst}; set BIOMARKER_INSTALL_DIR to the folder that holds env/ and ref/)")
        else ok('micromamba', "${inst}/env/bin/micromamba")
        def envs = ['fastp', 'star', 'salmon', 'samtools', 'subread', 'multiqc', 'pyqc', 'r_dge', 'r_classifier', 'r_pathway'] + (params.run_strand_check ? ['rseqc'] : [])
        def missingEnvs = envs.findAll { !new File("${inst}/env/prefix/${it}/bin").isDirectory() }
        if (missingEnvs) bad('conda environments', "missing in ${inst}/env/prefix: ${missingEnvs.join(', ')} (build them with scripts/build_envs.sh)")
        else ok('conda environments', "all ${envs.size()} tool environments found (${envs.join(', ')})")
        // ---- reference and STAR index
        def ref = params.ref_dir.toString()
        def relRef = ['GRCh38_primary_chr.fa', "gencode.${params.gencode_release}.primary_chr.annotation.gtf.gz", "gencode.${params.gencode_release}.transcripts.primary_chr.fa.gz",
                      'tx2gene.tsv', 'gene_id_to_symbol.tsv', 'salmon_index', 'star_index']
        def noRef = relRef.findAll { !new File("${ref}/${it}").exists() }
        if (noRef) bad('reference', "missing in ${ref}: ${noRef.join(', ')} (the reference is installed once, see docs/code_sources.md; do not point input_dir there)")
        else ok('reference', "genome, annotation, transcripts, tx2gene, symbols and Salmon index found in ${ref}")
        def si = ['Genome', 'SA', 'SAindex', 'genomeParameters.txt'].findAll { !new File("${ref}/star_index/${it}").isFile() || new File("${ref}/star_index/${it}").length() == 0 }
        if (si) bad('STAR index', "incomplete in ${ref}/star_index: ${si.join(', ')} missing or empty")
        else ok('STAR index', "Genome, SA, SAindex found (${gb(new File(ref + '/star_index/SA').length())} for SA)")
        // ---- locked model files
        def lk = new File("${params.code_dir}/lock/diagnostic_v1/sha256sums_final.txt")
        if (!lk.isFile()) bad('locked model', "${lk} not found: the code folder is incomplete (code_dir = ${params.code_dir})") else ok('locked model', 'lock/diagnostic_v1 found (hashes are verified when the run starts)')
        // ---- memory and CPUs
        def mem = null
        for (p in ['/sys/fs/cgroup/memory.max', '/sys/fs/cgroup/memory/memory.limit_in_bytes']) { def f = new File(p); if (f.exists()) { def t = f.text.trim(); mem = t.isLong() ? t.toLong() : null; break } }
        def minB = (params.min_mem_star_gb as long) * 1024L * 1024L * 1024L
        if (mem == null) bad('memory', 'the memory limit of this pod cannot be read (or is "max"): the aligner needs a known limit of at least ' + params.min_mem_star_gb + ' GiB')
        else if (mem < minB) bad('memory', "${gb(mem)} is below the ${params.min_mem_star_gb} GiB that the aligner (STAR) needs. Use a pod with 64 GB.")
        else if (mem < 60_000_000_000L) warn('memory', "${gb(mem)}: enough for the aligner, below the 64 GB recommended (tasks will queue)")
        else ok('memory', "${gb(mem)} (the aligner needs at least ${params.min_mem_star_gb} GiB)")
        def cpus = Runtime.runtime.availableProcessors()
        if (cpus < 16) warn('CPUs', "${cpus} CPUs; 16 are recommended (the run is slower, not wrong)") else ok('CPUs', "${cpus}")
    }

    def nBad = L.count { it[0] == '[MISSING]' }; def nWarn = L.count { it[0] == '[WARN]' }
    log.info "STARTUP CHECKS (mode ${mode}${workflow.stubRun && !params.check_only ? ', stub run: software and reference checks skipped' : ''})\n" +
        L.collect { "  ${it[0].padRight(9)} ${it[1].padRight(22)} ${it[2]}" }.join('\n') +
        "\n  => ${L.count { it[0] == '[ OK ]' }} ok, ${nWarn} warning(s), ${nBad} missing"
    if (nBad) error "STARTUP CHECKS FAILED: ${nBad} item(s) missing, nothing was run:\n" + L.findAll { it[0] == '[MISSING]' }.collect { "  - ${it[1]}: ${it[2]}" }.join('\n')
}
