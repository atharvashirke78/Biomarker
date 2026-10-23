#!/usr/bin/env nextflow
// Biomarker RNA-seq, diagnostic branch. One pipeline, two stages (docs/two_pipeline_design.md):
//   QUANT    FASTQ -> fastp -> (STAR -> samtools -> featureCounts) + Salmon/tximport -> per-sample QC -> cohort QC -> contract files + manifest
//   ANALYSIS counts matrix + sample sheet -> input checks -> DESeq2 -> pathways (ORA/GSEA)   [prognostic/predictive branches later: --analysis selector]
// nextflow run main.nf -entry QUANT | -entry ANALYSIS | (no -entry: QUANT then ANALYSIS, same as before the refactor)
// Per-sample steps never wait for other samples; cohort steps start after collect().
nextflow.enable.dsl = 2

include { MAKE_SAMPLE_MAP                                  } from '../modules/sample_map'
include { VERIFY_FASTQ                                   } from '../modules/verify_fastq'
include { SUBSAMPLE                                      } from '../modules/subsample'
include { FASTQC                                         } from '../modules/fastqc'
include { FASTP                                          } from '../modules/fastp'
include { REF_FETCH; REF_FILTER; REF_GENOME_FILTER; REF_TABLES; SALMON_INDEX } from '../modules/reference'
include { STAR_INDEX; STAR_ALIGN                         } from '../modules/star'
include { SAMTOOLS_SORT_INDEX; SAMTOOLS_STATS            } from '../modules/samtools'
include { GTF_TO_BED12; STRAND_CHECK; MERGE_STRAND_CHECK } from '../modules/strand_check'
include { SALMON_QUANT                                   } from '../modules/salmon'
include { FEATURECOUNTS                                  } from '../modules/featurecounts'
include { SAMPLE_QC_FLAGS; MERGE_SAMPLE_QC               } from '../modules/sample_qc'
include { MERGE_FEATURECOUNTS; SALMON_TXIMPORT           } from '../modules/merge_counts'
include { COHORT_QC                                      } from '../modules/cohort_qc'
include { DESEQ2 as DESEQ2_PRIMARY; DESEQ2 as DESEQ2_SAL; DE_CONCORDANCE } from '../modules/deseq2'
include { PATHWAY_ORA; PATHWAY_GSEA                      } from '../modules/pathway'
include { MAKE_CONTRACT; VALIDATE_ANALYSIS_INPUT as VALIDATE_PRIMARY; VALIDATE_ANALYSIS_INPUT as VALIDATE_SECONDARY; WRITE_MANIFEST as MANIFEST_QUANT; WRITE_MANIFEST as MANIFEST_ANALYSIS } from '../modules/contract'
include { MULTIQC; QC_COPY_FASTP                         } from '../modules/multiqc'
include { DISCOVER_INPUT_CHECK; DISCOVER_DESEQ2 as DISCOVER_DESEQ2_PRIMARY; DISCOVER_DESEQ2 as DISCOVER_DESEQ2_SECONDARY; CANDIDATE_RULES; FIT_MODEL; PACKAGE_MODEL } from '../modules/discover'
include { INTERNAL_VALIDATION } from '../modules/internal_validation'     // upgrade H3: additive internal validation, report-only (params.run_internal_validation)
include { folderInput; runStartInfo; startupChecks; FILTER_CLINICAL; CHECK_QC } from '../modules/folder_input'
include { LOCKED_SCORE                                  } from '../modules/lock_score'
include { REPORT                                        } from '../modules/report'
// check_vs_paper (Hosseini 2023) is validation-only and lives in validation/ (not part of QUANT or ANALYSIS)

// ---- QUANT sample sheet (csv). Required: sample,run,study,group. New columns: batch, baseline (TRUE/FALSE, one baseline group), count_source (ignored: QUANT fills it).
// Without batch/baseline the sheet is only accepted with --legacy_samplesheet true (batch = study, baseline = the group named 'normal').
// Group names are NOT fixed in code. Returns rows with batch and baseline resolved.
// QUANT sample sheet of this launch: the sheet built from --input_dir (-entry RUN, folderInput), else params.samplesheet. (params cannot be changed at run time, so the path is derived.)
def quantSheetPath() { params.input_dir ? "${params.output_dir}/00_input/samplesheet_quant.csv" : params.samplesheet }

def parseSamplesheet(String path) {
    def f = file(path)
    if (!f.exists()) error "Samplesheet not found: ${path}"
    def lines = f.readLines().findAll { it.trim() }
    def header = lines[0].split(',', -1)*.trim()
    def missing = ['sample', 'run', 'study', 'group'] - header
    if (missing) error "Samplesheet missing columns: ${missing.join(', ')}"
    if (('fastq_1' in header) != ('fastq_2' in header)) error "Samplesheet: fastq_1 and fastq_2 columns go together"
    def newcols = ['batch', 'baseline'] - header
    if (newcols && !params.legacy_samplesheet) error "Samplesheet lacks column(s) ${newcols.join(', ')}; add them or pass --legacy_samplesheet true (batch = study, baseline = group 'normal')"
    def rows = lines.drop(1).collect { l ->
        def c = l.split(',', -1)*.trim()
        if (c.size() != header.size()) error "Bad samplesheet row: ${l}"
        [header, c].transpose().collectEntries { k, v -> [(k): v] }
    }
    def dups = rows*.sample.countBy { it }.findAll { k, v -> v > 1 }.keySet()
    if (dups) error "Duplicate sample names: ${dups.join(', ')}"
    rows.each { r ->
        if (!r.group) error "Sample ${r.sample}: empty group"
        // folder entry (scripts/rnaseq.sh): FASTQ paths come from the sheet, IDs are free (checked by scripts/input_scan.py)
        if ('fastq_1' in header) { if (!r.fastq_1 || !r.fastq_2) error "Sample ${r.sample}: empty fastq_1/fastq_2" }
        else if (!(r.run ==~ /^[SED]RR[0-9]+$/)) error "Sample ${r.sample}: '${r.run}' is not an SRA run accession"
        if (!('batch' in header)) r.batch = r.study
        if (!('baseline' in header)) r.baseline = (r.group == 'normal') ? 'TRUE' : 'FALSE'
        r.baseline = r.baseline.toUpperCase()
        if (!(r.baseline in ['TRUE', 'FALSE'])) error "Sample ${r.sample}: baseline must be TRUE or FALSE, got '${r.baseline}'"
        if (!r.batch) error "Sample ${r.sample}: empty batch"
    }
    def bg = rows.findAll { it.baseline == 'TRUE' }*.group.unique()
    if (bg.size() != 1) error "Samplesheet: exactly one baseline group required, found ${bg ?: 'none'}" + (params.legacy_samplesheet && !('baseline' in header) ? " (legacy sheet: baseline = group named 'normal')" : "")
    rows.groupBy { it.group }.each { g, rs -> if (rs*.baseline.unique().size() > 1) error "Group ${g}: baseline flag differs between samples" }
    return rows
}

// [R1, R2] path strings of one sheet row: fastq_1/fastq_2 columns when present, else the legacy data/{normal,tumor}_sample/<run>_R{1,2}.fastq.gz
def fastqPair(Map r) {
    if (r.fastq_1) return [r.fastq_1, r.fastq_2]
    def d = r.baseline == 'TRUE' ? "${params.base_dir}/data/normal_sample" : "${params.base_dir}/data/tumor_sample"
    return ["${d}/${r.run}_R1.fastq.gz", "${d}/${r.run}_R2.fastq.gz"]
}

// cgroup memory limit in bytes (v2, then v1); null if unreadable or 'max'.
def cgroupMemBytes() {
    for (p in ['/sys/fs/cgroup/memory.max', '/sys/fs/cgroup/memory/memory.limit_in_bytes']) {
        def f = new File(p)
        if (f.exists()) { def t = f.text.trim(); return t.isLong() ? t.toLong() : null }
    }
    return null
}

// Memory guard (owner decision 2026-10-12): STAR route needs >= min_mem_star_gb. Below that, or limit unreadable, or run_star = false:
// ERROR, unless --allow_salmon_only true. The method actually used goes into the manifest and the report header.
def resolveMethod() {
    def mem = cgroupMemBytes()
    def minB = (params.min_mem_star_gb as long) * 1024L * 1024L * 1024L
    def why = null
    if (!params.run_star) why = "run_star = false"
    else if (mem == null) why = "the cgroup memory limit is unreadable (or 'max')"
    else if (mem < minB) why = "cgroup memory limit ${mem} bytes < ${params.min_mem_star_gb} GB (min_mem_star_gb)"
    if (why && !params.allow_salmon_only)
        error "STAR + featureCounts route cannot run: ${why}. Use a pod with >= ${params.min_mem_star_gb} GB, or pass --allow_salmon_only true to run the Salmon-only method (recorded in the manifest and report header)."
    def star = (why == null)
    return [mem: mem, mem_readable: mem != null, run_star: star, method: star ? 'STAR+featureCounts' : 'Salmon-only', route: star ? 'featurecounts' : 'salmon_tximport']
}

def paramsJson() {
    def m = params.collectEntries { k, v -> [(k): (v instanceof Number || v instanceof Boolean) ? v : v?.toString()] }
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(m))
}
def runMeta(String entry, Map extra) {
    def v = file("${projectDir}/../VERSION")
    def m = [pipeline_version: v.exists() ? v.text.trim() : 'unknown', entry: entry, run_name: workflow.runName, command_line: workflow.commandLine,
             start: workflow.start?.toString(), nextflow: workflow.nextflow.version.toString(), java: System.getProperty('java.runtime.version'),
             analysis: params.analysis, fastq_md5: params.manifest_fastq_md5] + extra
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(m))
}

// join helpers: tuples are [meta, ...]; join per sample id so samples never wait for each other
def keyL(ch) { ch.map { t -> [t[0].id] + t } }
def keyR(ch) { ch.map { t -> [t[0].id] + t.drop(1) } }
def unkey(ch) { ch.map { t -> t.drop(1) } }

workflow REFERENCE {
    main:
    def rel = params.gencode_release
    def tx  = file("${params.ref_dir}/gencode.${rel}.transcripts.primary_chr.fa.gz")
    def gt  = file("${params.ref_dir}/gencode.${rel}.primary_chr.annotation.gtf.gz")
    def gn  = file("${params.ref_dir}/GRCh38_primary_chr.fa")
    if (tx.exists() && gt.exists()) { ch_tx = Channel.value(tx); ch_gt = Channel.value(gt) }
    else { REF_FETCH(); REF_FILTER(REF_FETCH.out.gtf, REF_FETCH.out.transcripts); ch_tx = REF_FILTER.out.transcripts; ch_gt = REF_FILTER.out.gtf }
    if (gn.exists()) { ch_gn = Channel.value(gn) } else { REF_GENOME_FILTER(); ch_gn = REF_GENOME_FILTER.out.genome }

    def idx = file("${params.ref_dir}/salmon_index")
    if (idx.exists()) { ch_idx = Channel.value(idx) } else { SALMON_INDEX(ch_tx, ch_gn); ch_idx = SALMON_INDEX.out.index }
    def t2g = file("${params.ref_dir}/tx2gene.tsv"); def sym = file("${params.ref_dir}/gene_id_to_symbol.tsv")
    if (t2g.exists() && sym.exists()) { ch_t2g = Channel.value(t2g); ch_sym = Channel.value(sym) }
    else { REF_TABLES(ch_gt); ch_t2g = REF_TABLES.out.tx2gene; ch_sym = REF_TABLES.out.symbols }
    emit:
    salmon_index = ch_idx
    tx2gene = ch_t2g
    symbols = ch_sym
    genome = ch_gn
    gtf = ch_gt
}


// =====================================================================================================================
// QUANT stage: raw FASTQ -> contract files (<outdir>/contract/) + manifest_quant.json
// =====================================================================================================================
workflow QUANT_STAGE {
    main:
    def rm = resolveMethod()
    def run_star = rm.run_star
    log.info "cgroup memory limit: ${rm.mem} bytes (readable: ${rm.mem_readable}); method used: ${rm.method}"

    REFERENCE()
    def rows = parseSamplesheet(quantSheetPath())
    def noFile = file("${projectDir}/assets/NO_FILE")
    def noFileFc = file("${projectDir}/assets/NO_FILE_FC")
    // one resolved csv (sample,run,study,group,batch,baseline) is the single source for sample_map and the contract sheet
    def csv = (['sample,run,study,group,batch,baseline'] + rows.collect { r -> [r.sample, r.run, r.study, r.group, r.batch, r.baseline].join(',') }).join('\n') + '\n'
    ch_resolved = Channel.of(csv).collectFile(name: 'samplesheet_resolved.csv').first()
    MAKE_SAMPLE_MAP(ch_resolved)           // results/sample_map.tsv: srr, label, study, group (cohort outputs use the SRR as sample ID)
    smap = MAKE_SAMPLE_MAP.out.map

    // FASTQ paths: fastq_1/fastq_2 sheet columns (folder entry); else the legacy layout: baseline group in data/normal_sample, the other in data/tumor_sample
    ch_in = Channel.fromList(rows).map { r -> def p = fastqPair(r); tuple([id: r.sample, run: r.run, study: r.study, group: r.group], p[0], p[1]) }
    VERIFY_FASTQ(ch_in)            // per sample: released as soon as its own two files are complete
    ch_reads = VERIFY_FASTQ.out.reads
    if (params.subsample_pairs > 0) { SUBSAMPLE(ch_reads); ch_reads = SUBSAMPLE.out.reads }
    if (params.run_fastqc) { FASTQC(ch_reads) }

    FASTP(ch_reads)

    // ---------------- per-sample: Salmon
    SALMON_QUANT(FASTP.out.reads, REFERENCE.out.salmon_index)

    // ---------------- per-sample: STAR -> samtools -> featureCounts
    if (run_star) {
        def first_id = rows[0].sample
        def sidx = file("${params.ref_dir}/star_index")
        if (sidx.exists()) { ch_sidx = Channel.value(sidx) }
        else {
            ch_len_json = FASTP.out.json.filter { m, j -> m.id == first_id }.map { m, j -> j }
            STAR_INDEX(REFERENCE.out.genome, REFERENCE.out.gtf, ch_len_json)
            ch_sidx = STAR_INDEX.out.index.first()   // queue -> value channel: one index reused by all samples
        }
        STAR_ALIGN(FASTP.out.reads, ch_sidx)
        SAMTOOLS_SORT_INDEX(STAR_ALIGN.out.bam)
        SAMTOOLS_STATS(SAMTOOLS_SORT_INDEX.out.bam)
        if (params.run_strand_check) {        // flag-only second strandedness estimate (RSeQC); never changes counts, QC status or exclusions
            GTF_TO_BED12(REFERENCE.out.gtf)
            STRAND_CHECK(unkey(keyL(SAMTOOLS_SORT_INDEX.out.bam).join(keyR(SALMON_QUANT.out.quant))), GTF_TO_BED12.out.bed)
            MERGE_STRAND_CHECK(STRAND_CHECK.out.row.map { m, f -> f }.collect(), MAKE_SAMPLE_MAP.out.map)
            ch_strand = MERGE_STRAND_CHECK.out.table
        }
        FEATURECOUNTS(unkey(keyL(SAMTOOLS_SORT_INDEX.out.bam).join(keyR(SALMON_QUANT.out.quant))), REFERENCE.out.gtf)
        // SAMPLE_QC_FLAGS input: meta, star_log, fc_summary, quant, fastp_json
        ch_qc_in = keyL(SALMON_QUANT.out.quant)                        // [id, meta, quant]
            .join(keyR(STAR_ALIGN.out.log))                            // + star_log
            .join(keyR(FEATURECOUNTS.out.counts.map { m, c, s -> [m, s] }))   // + fc_summary
            .join(keyR(FASTP.out.json))                                // + fastp_json
            .map { id, m, q, sl, fs, fj -> tuple(m, sl, fs, q, fj) }
    } else {
        ch_qc_in = keyL(SALMON_QUANT.out.quant).join(keyR(FASTP.out.json))
            .map { id, m, q, fj -> tuple(m, noFile, noFileFc, q, fj) }   // two DIFFERENT placeholder names (same-name path inputs collide)
    }
    SAMPLE_QC_FLAGS(ch_qc_in)

    // ---------------- cohort steps: nothing below runs until all samples are done (collect)
    MERGE_SAMPLE_QC(SAMPLE_QC_FLAGS.out.row.map { m, f -> f }.collect(), smap)
    // status is read BY NAME (column 'status' of the merged sample_qc table, which has a header), never by position; the cohort steps wait for all samples anyway
    ch_ok = MERGE_SAMPLE_QC.out.table.splitCsv(header: true, sep: '\t').filter { r -> r.status != 'FAIL' }.map { r -> [r.srr, r.status] }
    ch_quant_ok = ch_ok.join(keyR(SALMON_QUANT.out.quant)).map { id, s, q -> q }.collect()
    SALMON_TXIMPORT(ch_quant_ok, REFERENCE.out.tx2gene, smap)

    if (run_star) {
        ch_fc_ok = ch_ok.join(keyR(FEATURECOUNTS.out.counts.map { m, c, s -> [m, c] })).map { id, s, c -> c }.collect()
        MERGE_FEATURECOUNTS(ch_fc_ok, smap)
        ch_primary_matrix = MERGE_FEATURECOUNTS.out.matrix
        COHORT_QC(MERGE_FEATURECOUNTS.out.matrix, smap, MERGE_SAMPLE_QC.out.table)
    } else {
        ch_primary_matrix = SALMON_TXIMPORT.out.counts
        COHORT_QC(SALMON_TXIMPORT.out.rds, smap, MERGE_SAMPLE_QC.out.table)
    }

    // ---------------- contract files + manifest
    MAKE_CONTRACT(ch_resolved, MERGE_SAMPLE_QC.out.table, COHORT_QC.out.flags, ch_primary_matrix,
                  SALMON_TXIMPORT.out.counts, SALMON_TXIMPORT.out.tpm, SALMON_TXIMPORT.out.rds, rm.route)
    def fq = rows.collect { r -> fastqPair(r) }.flatten()
    ch_mfiles = MAKE_CONTRACT.out.sheet.mix(MAKE_CONTRACT.out.qc, MAKE_CONTRACT.out.counts, ch_resolved).collect()
    MANIFEST_QUANT('quant', runMeta('QUANT', [method_used: rm.method, allow_salmon_only: params.allow_salmon_only, cgroup_memory_limit_bytes: rm.mem,
                                              cgroup_memory_limit_readable: rm.mem_readable, count_source: rm.route]),
                   paramsJson(), ch_mfiles, ['sample_sheet=sample_sheet.tsv', 'qc_table=qc_table.tsv', 'counts_matrix=counts_matrix.tsv', 'input_sample_sheet=samplesheet_resolved.csv'],
                   fq, noFile, true)

    // ---------------- MultiQC: once per results folder (that folder's inputs) and once over everything
    ch_fp  = FASTP.out.json_only.collect()
    ch_sal = SALMON_QUANT.out.quant.map { m, q -> q }.collect()
    ch_all = ch_fp.mix(ch_sal)
    mq = Channel.empty().mix(ch_fp.map { ['01_fastp', it] }, ch_sal.map { ['03_quant_salmon', it] })
    if (run_star) {
        ch_aln = STAR_ALIGN.out.log.map { m, f -> f }.mix(SAMTOOLS_STATS.out.stats.flatMap { m, a, b -> [a, b] }).collect()
        ch_fcs = FEATURECOUNTS.out.counts.map { m, c, s -> s }.collect()
        ch_all = ch_all.mix(ch_aln, ch_fcs)
        mq = mq.mix(ch_aln.map { ['02_alignment', it] }, ch_fcs.map { ['04_counts', it] })
    }
    ch_everything = ch_all.flatten().collect()
    mq = mq.mix(ch_everything.map { ['08_multiqc', it] })      // 01-04: one report per stage with its own files; 08: the combined report (no duplicate reports in 05-07)
    MULTIQC(mq, smap)
    QC_COPY_FASTP(FASTP.out.json_only.mix(FASTP.out.html).collect(), smap)

    emit:
    cohort_plots = COHORT_QC.out.pca.mix(COHORT_QC.out.heatmap, COHORT_QC.out.libsize, COHORT_QC.out.flags, COHORT_QC.out.summary, COHORT_QC.out.sex).mix(run_star && params.run_strand_check ? ch_strand : Channel.empty()).collect()   // plots + flags/summary/sex files (the report lists their WARNING lines)
    sample_qc = MERGE_SAMPLE_QC.out.table
    counts   = MAKE_CONTRACT.out.counts
    sal_counts = MAKE_CONTRACT.out.sal_counts
    rds      = MAKE_CONTRACT.out.rds
    sheet    = MAKE_CONTRACT.out.sheet
    qc       = MAKE_CONTRACT.out.qc
    symbols  = REFERENCE.out.symbols
    manifest = MANIFEST_QUANT.out.manifest
}

// =====================================================================================================================
// ANALYSIS stage: contract files -> input checks -> DESeq2 -> pathways. Never reads FASTQ, BAM or work/.
// route = featurecounts (counts_tsv, integers) | salmon_tximport (counts_rds)
// =====================================================================================================================
workflow ANALYSIS_STAGE {
    take:
    route
    counts_tsv
    counts_rds
    sheet
    symbols
    qc_table

    main:
    if (params.analysis != 'diagnostic') error "--analysis ${params.analysis}: only 'diagnostic' exists (prognostic/predictive branches are not implemented yet)"
    VALIDATE_PRIMARY('primary', counts_tsv, sheet, symbols, qc_table, route, 'false')
    DESEQ2_PRIMARY('deseq2', route.map { it == 'featurecounts' ? 'matrix' : 'tximport' }, counts_tsv, counts_rds, VALIDATE_PRIMARY.out.sheet, symbols)
    PATHWAY_ORA(DESEQ2_PRIMARY.out.results)
    PATHWAY_GSEA(DESEQ2_PRIMARY.out.results)

    emit:
    results = DESEQ2_PRIMARY.out.results
    sheet   = VALIDATE_PRIMARY.out.sheet
    ora     = PATHWAY_ORA.out.tables
    gsea    = PATHWAY_GSEA.out.tables
    dge_plots     = DESEQ2_PRIMARY.out.volcano.mix(DESEQ2_PRIMARY.out.heatmap, DESEQ2_PRIMARY.out.ma).collect()
    pathway_plots = PATHWAY_ORA.out.plots.mix(PATHWAY_GSEA.out.plots).flatten().ifEmpty { file("${projectDir}/assets/NO_FILE") }.collect()
}

// optional cross-check route (default run, STAR route only): Salmon/tximport DESeq2 + concordance with the primary route
workflow SECONDARY_ROUTE {
    take:
    primary_results
    sal_counts_tsv
    sal_rds
    sheet
    symbols
    qc_table

    main:
    VALIDATE_SECONDARY('secondary', sal_counts_tsv, sheet, symbols, qc_table, 'salmon_tximport', 'true')
    DESEQ2_SAL('deseq2_salmon', 'tximport', sal_counts_tsv, sal_rds, VALIDATE_SECONDARY.out.sheet, symbols)
    DE_CONCORDANCE(primary_results, DESEQ2_SAL.out.results)
}

// =====================================================================================================================
// Entry points
// =====================================================================================================================
workflow QUANT {
    QUANT_STAGE()
}

// -entry ANALYSIS --counts_matrix <counts_matrix.tsv> --sample_sheet <sample_sheet.tsv> [--count_source featurecounts|salmon_tximport] [--counts_rds <rds>] [--qc_table <qc_table.tsv>] [--quant_manifest <manifest_quant.json>]
workflow ANALYSIS {
    main:
    def need = ['counts_matrix', 'sample_sheet']
    need.each { if (!params[it]) error "-entry ANALYSIS needs --${it}" }
    def route = params.count_source
    if (!(route in ['featurecounts', 'salmon_tximport'])) error "--count_source must be featurecounts or salmon_tximport, got '${route}'"
    def noFile = file("${projectDir}/assets/NO_FILE")
    def tsv = file(params.counts_matrix); def sheet = file(params.sample_sheet)
    if (!tsv.exists()) error "counts matrix not found: ${params.counts_matrix}"
    if (!sheet.exists()) error "sample sheet not found: ${params.sample_sheet}"
    def rds = params.counts_rds ? file(params.counts_rds) : noFile
    if (route == 'salmon_tximport' && !rds.exists()) error "count_source salmon_tximport needs --counts_rds (the tximport object; --counts_matrix is then the matching counts TSV)"
    def qc = params.qc_table ? file(params.qc_table) : noFile
    def qm = params.quant_manifest ? file(params.quant_manifest) : noFile
    def sym = file(params.symbols ?: "${params.ref_dir}/gene_id_to_symbol.tsv")
    if (!sym.exists()) error "gene_id_to_symbol table not found: ${sym} (pass --symbols)"

    ANALYSIS_STAGE(Channel.value(route), Channel.value(tsv), Channel.value(rds), Channel.value(sheet), Channel.value(sym), Channel.value(qc))
    ch_files = ANALYSIS_STAGE.out.sheet.mix(ANALYSIS_STAGE.out.results).collect()
    MANIFEST_ANALYSIS('analysis', runMeta('ANALYSIS', [count_source: route, method_used: 'n/a (ANALYSIS only; see quant_manifest if given)']), paramsJson(),
                      Channel.value(tsv).mix(Channel.value(sheet), ch_files.flatten()).collect(),
                      ['counts_matrix=' + tsv.name, 'sample_sheet=' + sheet.name, 'validated_sheet=sample_sheet.validated.tsv', 'results_all=deseq2_results_all.tsv'],
                      [], qm, false)
}

// -entry DISCOVER --question diagnostic --counts_matrix <featureCounts counts.tsv> --clinical <clinical.tsv> [--holdout <ids.txt>]
//                 [--counts_matrix_secondary <salmon counts.tsv> --counts_rds_secondary <tximport.rds>]   (Phase 1, 2026-10-15; plan: results_rc1prep_19/pipeline_plan_20261015/PLAN.md)
// Input contract (stops on bad input) -> DESeq2 per route -> candidate rules 1-3b -> diagnostic model fit -> packaging. Never writes into lock/.
workflow DISCOVER {
    main:
    if (!params.question) error "-entry DISCOVER needs --question diagnostic|prognostic|predictive"
    if (params.question != 'diagnostic') error "--question ${params.question}: only 'diagnostic' is implemented (prognostic/predictive: Phase 2)"
    ['counts_matrix', 'clinical'].each { if (!params[it]) error "-entry DISCOVER needs --${it}" }
    def dual = params.counts_matrix_secondary as boolean
    if (dual && !params.counts_rds_secondary) error "--counts_matrix_secondary needs --counts_rds_secondary (the matching tximport .rds)"
    def noFile = file("${projectDir}/assets/NO_FILE")
    def counts = file(params.counts_matrix); def clinical = file(params.clinical)
    def sym = file(params.symbols ?: "${params.ref_dir}/gene_id_to_symbol.tsv")
    [counts, clinical, sym].each { if (!it.exists()) error "input not found: ${it}" }
    def sec = dual ? file(params.counts_matrix_secondary).toString() : ''
    def ho = params.holdout ? file(params.holdout).toString() : ''

    DISCOVER_INPUT_CHECK(counts, clinical, sym, sec, ho)
    sheet = DISCOVER_INPUT_CHECK.out.sheet
    ctr = DISCOVER_INPUT_CHECK.out.counts_training
    DISCOVER_DESEQ2_PRIMARY('deseq2', 'matrix', ctr, noFile, sheet, sym)
    if (dual) {
        DISCOVER_DESEQ2_SECONDARY('deseq2_salmon', 'tximport', file(params.counts_matrix_secondary), file(params.counts_rds_secondary), sheet, sym)
        de2 = DISCOVER_DESEQ2_SECONDARY.out.results
    } else {
        de2 = Channel.value(noFile)
    }
    def rules = file("${params.code_dir}/config/candidate_rules.yaml"); def thr = file("${params.code_dir}/conf/thresholds.config"); def qcfg = file(params.questions_config)
    CANDIDATE_RULES(ctr, sheet, DISCOVER_DESEQ2_PRIMARY.out.results, de2, rules, thr, qcfg)
    FIT_MODEL(ctr, sheet, CANDIDATE_RULES.out.table, DISCOVER_DESEQ2_PRIMARY.out.results, qcfg, CANDIDATE_RULES.out.note)
    PACKAGE_MODEL(FIT_MODEL.out.genes, FIT_MODEL.out.spec, FIT_MODEL.out.folds, FIT_MODEL.out.scores, FIT_MODEL.out.perf,
                  DISCOVER_INPUT_CHECK.out.checks, CANDIDATE_RULES.out.note, thr, qcfg, counts, clinical, CANDIDATE_RULES.out.table)
}

// =====================================================================================================================
// -entry RUN (and the default when --input_dir is given): ONE workflow, FASTQ folder in, biomarker list and report out.
//   nextflow run nextflow/main.nf --input_dir <folder> [--output_dir <folder>] [--mode discover|apply] [--cohort development|validation] [-resume]
//   discover (default): FOLDER_INPUT -> QUANT_STAGE -> ANALYSIS_STAGE (DESeq2 + pathways) -> DISCOVER (candidate rules, LASSO model) -> REPORT
//   apply:              FOLDER_INPUT -> QUANT_STAGE -> LOCKED_SCORE (lock/diagnostic_v1) -> REPORT
// Everything is written under output_dir (default <input_dir>/Result): numbered folders per stage, 10_report/, logs/, work/, run_info.json, params.json.
// =====================================================================================================================
workflow RUN_PIPELINE {
    main:
    def mode = params.mode
    def sa = params.stop_after
    if (sa && !(sa in ['quant', 'analysis', 'discover'])) error "--stop_after ${sa}: allowed values are quant, analysis, discover"
    if (sa in ['analysis', 'discover'] && mode != 'discover')
        error "--stop_after ${sa} applies to --mode discover only (apply mode has two parts, quantification and scoring: use --stop_after quant to stop before the locked-model scoring)"
    if (params.analysis != 'diagnostic') error "--analysis ${params.analysis}: only 'diagnostic' exists (prognostic/predictive branches are not implemented yet)"
    if (params.question != 'diagnostic') error "--question ${params.question}: only 'diagnostic' is implemented (prognostic/predictive: Phase 2)"
    if (mode == 'apply' && !(params.cohort in ['development', 'validation']))
        error "--mode apply needs --cohort development|validation (written into every score output; development = optimistic, not validation)"
    def fi = folderInput(mode)                       // input checks: stops here, before any compute, with a clear message
    def rm = resolveMethod()                         // memory guard (never bypassed here)
    if (!rm.run_star) error "-entry RUN needs the STAR + featureCounts route (${rm.method} is not supported by the one-workflow entry)"
    def info = runStartInfo(mode, paramsJson())      // lock hash check (stops on mismatch), run_info.json, params.json
    def noFile = file("${projectDir}/assets/NO_FILE")

    QUANT_STAGE()
    def q = QUANT_STAGE.out
    CHECK_QC(q.sample_qc, mode == 'discover' ? params.analysis_min_samples_per_group : 0)

    def scanNames = ['scan_files.tsv', 'problem_files.tsv', 'metadata.used.csv', 'metadata.resolved.tsv', 'group_mapping.tsv', 'samplesheet_quant.csv',
                     'clinical_discover.tsv', 'sample_sheet_apply.tsv', 'input_summary.json']
    ch_scan = Channel.fromList(scanNames.collect { file("${fi.dir}/${it}") }.findAll { it.exists() }).mix(CHECK_QC.out.excluded)
    def doneInfo = info + [stop_after: sa, attempts: info.attempts.take(info.attempts.size() - 1) + [info.attempts.last() + [status: 'completed', exit: 0]]]
    ch_info = Channel.of(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(doneInfo))).collectFile(name: 'run_info.json').first()

    // report inputs; stages that are not run (stop_after) arrive as the NO_FILE placeholder and are left out of the report
    ch_input = ch_scan.collect()
    ch_de = Channel.value(noFile); ch_cand = Channel.value(noFile); ch_model = Channel.value(noFile); ch_scores = Channel.value(noFile)
    ch_dge_plots = Channel.value(noFile); ch_path_plots = Channel.value(noFile); ch_valid = Channel.value(noFile)

    if (mode == 'discover') {
        if (sa != 'quant') {
            ANALYSIS_STAGE(Channel.value(rm.route), q.counts, q.rds, q.sheet, q.symbols, q.qc)
            ch_dge_plots = ANALYSIS_STAGE.out.dge_plots
            ch_path_plots = ANALYSIS_STAGE.out.pathway_plots
        }
        if (!sa || sa == 'discover') {          // the discovery chain runs unless the run is stopped after quant or analysis (stop_after discover = the whole discover mode)
            FILTER_CLINICAL(Channel.value(fi.clinical), q.counts)
            DISCOVER_INPUT_CHECK(q.counts, FILTER_CLINICAL.out.clinical, q.symbols, q.sal_counts.map { it.toString() }, '')
            DISCOVER_DESEQ2_PRIMARY('deseq2', 'matrix', DISCOVER_INPUT_CHECK.out.counts_training, noFile, DISCOVER_INPUT_CHECK.out.sheet, q.symbols)
            DISCOVER_DESEQ2_SECONDARY('deseq2_salmon', 'tximport', q.sal_counts, q.rds, DISCOVER_INPUT_CHECK.out.sheet, q.symbols)
            def rules = file("${params.code_dir}/config/candidate_rules.yaml"); def thr = file("${params.code_dir}/conf/thresholds.config"); def qcfg = file(params.questions_config)
            CANDIDATE_RULES(DISCOVER_INPUT_CHECK.out.counts_training, DISCOVER_INPUT_CHECK.out.sheet, DISCOVER_DESEQ2_PRIMARY.out.results,
                            DISCOVER_DESEQ2_SECONDARY.out.results, rules, thr, qcfg)
            FIT_MODEL(DISCOVER_INPUT_CHECK.out.counts_training, DISCOVER_INPUT_CHECK.out.sheet, CANDIDATE_RULES.out.table, DISCOVER_DESEQ2_PRIMARY.out.results,
                      qcfg, CANDIDATE_RULES.out.note)
            PACKAGE_MODEL(FIT_MODEL.out.genes, FIT_MODEL.out.spec, FIT_MODEL.out.folds, FIT_MODEL.out.scores, FIT_MODEL.out.perf,
                          DISCOVER_INPUT_CHECK.out.checks, CANDIDATE_RULES.out.note, thr, qcfg, q.counts, FILTER_CLINICAL.out.clinical, CANDIDATE_RULES.out.table)
            ch_input = ch_scan.mix(FILTER_CLINICAL.out.excluded).collect()
            ch_de    = DISCOVER_DESEQ2_PRIMARY.out.summary.mix(DISCOVER_DESEQ2_PRIMARY.out[3].flatten()).collect()     // [3] = the volcano/MA/heatmap PNGs
            ch_cand  = CANDIDATE_RULES.out.table.mix(CANDIDATE_RULES.out[2], CANDIDATE_RULES.out.note, CANDIDATE_RULES.out.refits).collect()   // [2] = funnel_counts.tsv
            ch_model = FIT_MODEL.out.genes.mix(FIT_MODEL.out.spec, FIT_MODEL.out.folds, FIT_MODEL.out.scores, FIT_MODEL.out.perf, FIT_MODEL.out.sens_spec,
                                               PACKAGE_MODEL.out[0], PACKAGE_MODEL.out[1]).collect()                  // report waits for the packaged model
            // ---- upgrade H3 (additive, report-only): internal validation of the finished model -> <outdir>/09_discover/04_validation. No-op when run_internal_validation = false.
            if (params.run_internal_validation) {
                INTERNAL_VALIDATION(DISCOVER_INPUT_CHECK.out.counts_training, DISCOVER_INPUT_CHECK.out.sheet, CANDIDATE_RULES.out.table,
                                    FIT_MODEL.out.genes, FIT_MODEL.out.spec, FIT_MODEL.out.scores)
                ch_valid = INTERNAL_VALIDATION.out.summary.mix(INTERNAL_VALIDATION.out.freq, INTERNAL_VALIDATION.out.md).collect()   // the report waits for the step and shows its files
            }
            // ---- end upgrade H3
        }
    } else if (sa != 'quant') {
        def evaluate = fi.summary.unknown_samples.isEmpty() && fi.summary.groups.size() == 2
        LOCKED_SCORE(q.counts, Channel.value(fi.apply_sheet), q.qc, Channel.value(file("${params.code_dir}/lock/diagnostic_v1")), params.cohort, evaluate)
        ch_scores = LOCKED_SCORE.out.scores
    }
    REPORT(ch_info, ch_input, q.qc, ch_de, ch_cand, ch_model, ch_scores, q.cohort_plots, ch_dge_plots, ch_path_plots, q.manifest, ch_valid)
}

// The entry point: startup checks (plain-language list of what is missing), then the pipeline. --check_only true prints the checks and stops.
workflow RUN {
    main:
    startupChecks(params.mode)
    if (!params.check_only) { RUN_PIPELINE() }
}

// default (no -entry): with --input_dir the one-workflow entry RUN; with --legacy_default true the old run (QUANT, then ANALYSIS fed by the contract files, unchanged); otherwise a clear stop
workflow {
    if (params.input_dir) {
        RUN()
    } else if (!params.legacy_default) {
        error "NO INPUT FOLDER SET. Nothing was run. Set input_dir in the USER SETTINGS block at the top of nextflow/nextflow.config, or start with --input_dir <folder with the paired FASTQ files and metadata.csv>. The pipeline never uses a default folder or test data."
    } else {
        def rm = resolveMethod()
        QUANT_STAGE()
        ANALYSIS_STAGE(Channel.value(rm.route), QUANT_STAGE.out.counts, QUANT_STAGE.out.rds, QUANT_STAGE.out.sheet, QUANT_STAGE.out.symbols, QUANT_STAGE.out.qc)
        if (rm.run_star) {
            SECONDARY_ROUTE(ANALYSIS_STAGE.out.results, QUANT_STAGE.out.sal_counts, QUANT_STAGE.out.rds,
                            QUANT_STAGE.out.sheet, QUANT_STAGE.out.symbols, QUANT_STAGE.out.qc)   // original contract sheet (the validated copy of the primary route is not reused)
        }
        ch_files = QUANT_STAGE.out.sheet.mix(ANALYSIS_STAGE.out.sheet, ANALYSIS_STAGE.out.results, QUANT_STAGE.out.counts).collect()
        MANIFEST_ANALYSIS('analysis', runMeta('QUANT+ANALYSIS', [count_source: rm.route, method_used: rm.method, allow_salmon_only: params.allow_salmon_only, cgroup_memory_limit_bytes: rm.mem]), paramsJson(),
                          ch_files, ['contract_sample_sheet=sample_sheet.tsv', 'validated_sheet=sample_sheet.validated.tsv', 'results_all=deseq2_results_all.tsv', 'counts_matrix=counts_matrix.tsv'],
                          [], QUANT_STAGE.out.manifest, false)
    }
}

workflow.onComplete {
    def em = (workflow.commandLine =~ /-entry\s+(\w+)/)
    def entry = em.find() ? em.group(1) : 'DEFAULT'
    // a run that stopped at the input check must leave nothing behind: no results folder in the checkout, none under a mistyped input path
    if (entry == 'DEFAULT' && !params.input_dir && !params.legacy_default) return
    if ((entry == 'RUN' || (entry == 'DEFAULT' && params.input_dir)) && (!params.input_dir || !params.output_dir || !new File(params.input_dir.toString()).isDirectory())) return
    def method = 'n/a'
    try { method = resolveMethod().method } catch (Exception e) { method = 'n/a (not resolved in this entry/pod)' }
    if (entry == 'RUN' || (entry == 'DEFAULT' && params.input_dir)) {
        def od = new File(params.output_dir.toString()); od.mkdirs()
        def ok = workflow.success
        def lockEnd = 'not checked'
        try {
            def pr = new ProcessBuilder(['sha256sum', '-c', 'lock/diagnostic_v1/sha256sums_final.txt']).directory(new File(params.code_dir.toString())).redirectErrorStream(true).start()
            def lt = pr.inputStream.text; lockEnd = (pr.waitFor() == 0 ? 'OK' : 'MISMATCH') + " (${lt.readLines().count { it.endsWith(': OK') }} of ${lt.readLines().findAll { it.contains(': ') }.size()})"
        } catch (Exception e) { }
        def jf = new File(od, 'run_info.json')
        if (jf.exists()) {
            def info = new groovy.json.JsonSlurper().parse(jf)
            def last = info.attempts.last()
            last.status = ok ? 'completed' : 'failed'; last.exit = workflow.exitStatus; last.end = workflow.complete?.toString(); last.lock_check_end = lockEnd
            def mqf = new File(od, 'contract/manifest_quant.json')       // reference files: path, size, md5 (STAR index: md5 of genomeParameters.txt + sizes) copied from the quantification manifest; no new hashing
            if (mqf.exists()) {
                def mq = new groovy.json.JsonSlurper().parse(mqf)
                info.reference = [gencode_release: mq.parameters?.gencode_release ?: params.gencode_release, files: mq.references,
                                  note: 'copied from contract/manifest_quant.json: md5 and sizes (md5 only; the whole STAR/Salmon index is not hashed)']
            }
            jf.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
        }
        new File(od, 'RESULTS_INDEX_RUN.md').text = """# Results index (run: ${workflow.runName}, entry: RUN, mode: ${params.mode}, ${ok ? 'completed' : 'FAILED'}, ${workflow.complete})
Case vs reference groups as given in the metadata. Development / technical output only; not clinically validated. Lock check at the end of the run: ${lockEnd}.
- 00_input/        folder scan: scan_files.tsv, problem_files.tsv, metadata.resolved.tsv, group_mapping.tsv, samplesheet_quant.csv, clinical tables, input_summary.json (notes, e.g. batch default)
- contract/        QUANT -> ANALYSIS files (counts_matrix.tsv, qc_table.tsv, sample_sheet.tsv, ...)
- 01_fastp/ 02_alignment/ 03_quant_salmon/ 04_counts/ 05_cohort_qc/ 08_multiqc/   trimming, STAR + samtools, Salmon, featureCounts, QC, MultiQC (01-04 each have the MultiQC of their own stage in qc/; 08_multiqc has the combined report)
- 06_dge/ 07_pathway/   DESeq2 on all QC-pass samples and ORA/GSEA (discover mode)
- 09_discover/     discover mode: 00_input (checks), 01_de, 02_candidates (candidate_table.tsv = biomarker candidates), 03_model (LASSO model, MODEL_MANIFEST.tsv, DISCOVER_SUMMARY.md)
- 09_locked_model/ apply mode: scores/ (scores.tsv, evaluation.tsv when labels are complete)
- 10_report/       report/report.html + report.md, tables/, plots/, manifest.json, tool_versions.tsv
- logs/            Nextflow trace, report and timeline;  run_info.json, params.json, tool_versions.tsv: run record;  work/: Nextflow work folder (needed for -resume)
"""
        return
    }
    if (entry == 'DISCOVER') {
        def di = new File("${params.outdir}/RESULTS_INDEX_DISCOVER.md"); di.parentFile.mkdirs()
        di.text = """# Results index (run: ${workflow.runName}, entry: DISCOVER, question: ${params.question}, ${workflow.success ? 'completed' : 'FAILED'}, ${workflow.complete})
Case vs reference groups as given in the metadata. Development / technical output only; not clinically validated. KEGG: ${params.run_kegg ? 'on' : 'off (default; licence unresolved)'}.
- 00_input/      input_checks.json (PASS/WARN/FAIL per check), sample_sheet.discover.tsv, training/hold-out sample lists
- 01_de/         DESeq2 per route (deseq2_*, deseq2_salmon_* when a second route is given)
- 02_candidates/ candidate rules 1-3b: candidate_table.tsv, funnel_counts.tsv, candidate_rules_method_note.txt (routes, skipped rules)
- 03_model/      model_genes.tsv, model_spec.txt, inner_fold_cutoffs.tsv, development_scores/performance (OPTIMISTIC), MODEL_MANIFEST.tsv, DISCOVER_SUMMARY.md
- */*.provenance.txt  command, tool versions and sha256 of inputs and outputs per process
"""
        return
    }
    def idx = new File("${params.outdir}/RESULTS_INDEX${entry == 'ANALYSIS' ? '_ANALYSIS' : ''}.md")
    idx.parentFile.mkdirs()
    idx.text = """# Results index (run: ${workflow.runName}, entry: ${entry}, ${workflow.success ? 'completed' : 'FAILED'}, ${workflow.complete})
Method used: ${method}   |   pipeline version: ${file("${projectDir}/../VERSION").text.trim()}   |   KEGG licence for commercial use is unresolved (see docs)
- contract/                 QUANT -> ANALYSIS files: counts_matrix.tsv, counts_salmon_tximport.tsv, salmon_tpm.tsv, salmon_tximport.rds, sample_sheet.tsv, qc_table.tsv (with reason), manifest_quant.json, manifest_analysis.json
- 00_analysis_input/        ANALYSIS input checks: <primary|secondary>_analysis_input_check.json and the validated sample sheet
- sample_map.tsv            srr, label, study, group (all per-sample folders/files and cohort sample IDs use the SRR run accession)
- 01_fastp/<SRR>/           fastp.html + fastp.json per sample (qc/: all samples' fastp files + MultiQC of fastp)
- 02_alignment/<SRR>/       <SRR>.sorted.bam + index, STAR Log.final.out, samtools flagstat and stats (BAMs only here)
- 03_quant_salmon/<SRR>/    Salmon output (quant.sf, aux_info/meta_info.json, lib_format_counts.json)
- 04_counts/                <SRR>/<SRR>.featureCounts.txt per sample, counts_featurecounts.tsv, counts_salmon_tximport.tsv, salmon_tpm.tsv, sample_qc.tsv (srr, label, study, group + PASS/WARN/FAIL + reason)
- 05_cohort_qc/             PCA (group colour, study shape), sample correlation heatmap, library sizes, MAD flags (flag only)
- 06_dge/                   DESeq2 (~ batch + group, reference = baseline group of the sample sheet): all/up/down tables, volcano, MA, top-50 heatmap; Salmon route (deseq2_salmon_*); concordance
- 07_pathway/               ORA and GSEA
- 08_multiqc/               combined MultiQC report over all inputs
- <folder>/qc/fastp/<SRR>.fastp.{html,json} fastp files of all samples; 01-04: <folder>/qc/multiqc_report.html = MultiQC of that stage only; 08_multiqc/multiqc_report.html = combined report over all stages (05-07 have no MultiQC report of their own)
"""
}
