#!/usr/bin/env python3
"""Report folder of a scripts/rnaseq.sh run (plan section 7.4). Reads the run folder only; never writes into lock/ or the input folder.

usage (project root): env/prefix/pyqc/bin/python3 -I scripts/make_report.py --run <run folder>
Writes into the run folder (each item only if it does not exist yet; never overwritten):
  qc/ align/ counts/ sample_qc/ deg/ candidates/ scores/   relative links to the pipeline output folders that exist
  model/diagnostic_v2_candidate_<run id>/                   copy of the DISCOVER model files (a CANDIDATE, never written into lock/)
  tables/ plots/                                            copies of the key tables and the plots
  tool_versions.tsv                                         copy of docs/tool_versions.tsv (the environments the pipeline uses)
  manifest.json                                             inputs, config, lock files and every report/table/model file with sha256
  report/report.md + report/report.html                     plain-language summary
report/ is written only for a completed attempt; a failed attempt writes report_partial_<UTC time>/; a later completed attempt when report/
exists writes report_<UTC time>/; manifest_<UTC time>.json likewise. Nothing is replaced."""
import argparse
import csv
import datetime
import glob
import hashlib
import html
import json
import os
import shutil
import subprocess

ap = argparse.ArgumentParser()
ap.add_argument('--run', required=True)
a = ap.parse_args()
R = os.path.realpath(a.run)
P = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
J = lambda *p: os.path.join(R, *p)
TS = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S')
info = json.load(open(J('run_info.json')))
RID, MODE = info['run_id'], info['mode']


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def rows(p):
    with open(p, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def md_table(rs, cols=None, limit=40):
    if not rs:
        return '_(none)_\n'
    cols = cols or list(rs[0])
    out = ['| ' + ' | '.join(cols) + ' |', '|' + '--|' * len(cols)]
    out += ['| ' + ' | '.join(str(r.get(c, '')) for c in cols) + ' |' for r in rs[:limit]]
    if len(rs) > limit:
        out.append(f'\n_{len(rs) - limit} more rows in the table file_')
    return '\n'.join(out) + '\n'


def stage_exit(name):
    p = J(f'nf_{name}', 'exit_code.txt')
    return open(p).read().strip() if os.path.exists(p) else None


# ---- 1. folder layout (links only where the target exists)
links = {'qc': 'pipeline/quant/08_multiqc', 'align': 'pipeline/quant/02_alignment', 'counts': 'pipeline/quant/04_counts',
         'sample_qc': 'pipeline/quant/05_cohort_qc', 'deg': 'pipeline/discover/01_de', 'candidates': 'pipeline/discover/02_candidates',
         'scores': 'pipeline/apply/scores'}
made = []
for k, t in links.items():
    if os.path.isdir(J(t)) and not os.path.lexists(J(k)):
        os.symlink(t, J(k))
        made.append(k)
mdir = J('model', f'diagnostic_v2_candidate_{RID}')
src_model = J('pipeline', 'discover', '03_model')
if os.path.isdir(src_model) and not os.path.exists(mdir):
    shutil.copytree(src_model, mdir)
    open(os.path.join(mdir, 'CANDIDATE_NOT_LOCKED.txt'), 'w').write(
        f'Candidate model of run {RID}. NOT a locked model: it is not in lock/ and has no tag. Locking is a separate owner step.\n')
tab_src = {'group_mapping.tsv': 'input/group_mapping.tsv', 'metadata.resolved.tsv': 'input/metadata.resolved.tsv',
           'problem_files.tsv': 'input/problem_files.tsv', 'sample_qc.tsv': 'pipeline/quant/04_counts/sample_qc.tsv',
           'excluded_after_quant.tsv': 'input/excluded_after_quant.tsv',
           'funnel_counts.tsv': 'pipeline/discover/02_candidates/funnel_counts.tsv', 'candidate_table.tsv': 'pipeline/discover/02_candidates/candidate_table.tsv',
           'model_genes.tsv': 'pipeline/discover/03_model/model_genes.tsv', 'development_performance.tsv': 'pipeline/discover/03_model/development_performance.tsv',
           'deseq2_summary.json': 'pipeline/discover/01_de/deseq2_summary.json', 'scores.tsv': 'pipeline/apply/scores/scores.tsv',
           'evaluation.tsv': 'pipeline/apply/scores/evaluation.tsv', 'evaluation_sens_spec.tsv': 'pipeline/apply/scores/evaluation_sens_spec.tsv', 'development_sens_spec.tsv': 'pipeline/discover/03_model/development_sens_spec.tsv', 'excluded_samples.tsv': 'pipeline/apply/scores/excluded_samples.tsv'}
for d in ('tables', 'plots'):
    os.makedirs(J(d), exist_ok=True)
for dst, src in tab_src.items():
    if os.path.exists(J(src)) and not os.path.exists(J('tables', dst)):
        shutil.copyfile(J(src), J('tables', dst))
plot_src = {}        # copied plot name -> stage folder it came from (used for the Plots section of the report)
for png in sorted(glob.glob(J('pipeline', '*', '*', '*.png'))):
    stage = os.path.relpath(os.path.dirname(png), J('pipeline'))
    name = os.path.basename(png)
    if os.path.exists(J('plots', name)):                      # same file name from another stage: keep both, prefix the stage folder
        name = os.path.basename(stage) + '_' + name
    if not os.path.exists(J('plots', name)):
        shutil.copyfile(png, J('plots', name))
    plot_src[name] = stage
if not os.path.exists(J('tool_versions.tsv')):
    shutil.copyfile(os.path.join(P, 'docs', 'tool_versions.tsv'), J('tool_versions.tsv'))

# ---- 2. report text
inp = json.load(open(J('input', 'input_summary.json'))) if os.path.exists(J('input', 'input_summary.json')) else None
att = info['attempts']
cohort = None                      # apply mode: the cohort label the locked-model scoring wrote into every score row
_sc = J('pipeline', 'apply', 'scores', 'scores.tsv')
if os.path.exists(_sc):
    _r = rows(_sc)
    cohort = _r[0].get('cohort') if _r else None
if cohort == 'development':
    banner = ('> **COHORT: DEVELOPMENT.** These samples were used to build the locked model: this result is NOT the independent validation and it is optimistic. '
              'Development / technical output only; not clinically validated.')
elif cohort == 'validation':
    banner = '> **COHORT: VALIDATION** (independent samples; locked model applied once). Technical output; not clinically validated.'
else:
    banner = '> Development / technical output only; not clinically validated. Development performance is OPTIMISTIC (genes and cutoff chosen on the same samples).'
_msp = J('pipeline', 'discover', '03_model', 'model_spec.txt')
MODEL_FITTED = os.path.exists(_msp) and 'NOT FITTED' not in open(_msp).read()      # gene selection (R4 banner) happened only when a model was fitted
PLAIN = '> Gene selection was done on all development samples; development performance is optimistic; only independent validation counts.'
VALID = '> These samples were not used to build the locked model. The locked model (lock/diagnostic_v1) was applied without any change. Quality-control exclusions are listed below.'
L = [f'# Run report: {RID}', '',
     f'**Mode:** {MODE} | **input:** `{info["input_dir"]}`' + (f' | **counts:** `{info["counts"]}`' if info.get('counts') else ' | **data:** FASTQ'),
     f'**Code:** commit `{info["git_commit"][:7]}` ({info["git_branch"]}), {len(info["uncommitted_files"])} uncommitted file(s) at start | '
     f'**lock/diagnostic_v1:** {info["lock_diagnostic_v1"]} | **pod:** {info["pod"]["nproc"]} CPUs, cgroup memory {info["pod"]["cgroup_memory_max"]} B',
     f'**Attempts:** ' + '; '.join(f'{x["kind"]} {x["start"]}' + (f' -> {x.get("status")} (exit {x.get("exit")})' if 'status' in x else ' (running)') for x in att),
     '', '**Intended use:** case-vs-reference biomarker discovery from RNA-seq (groups as given in metadata.csv); research and development use only; not clinically validated.', '', banner] + (
         ['', VALID] if cohort == 'validation' else ['', PLAIN] if ((MODE == 'discover' and MODEL_FITTED) or cohort == 'development') else [])
inp0 = json.load(open(J('input', 'input_summary.json'))) if os.path.exists(J('input', 'input_summary.json')) else None
SMALL_COHORT = 10      # PROVISIONAL: discover mode with a group smaller than this is a mechanics test
if MODE == 'discover' and inp0 and inp0.get('groups') and min(inp0['groups'].values()) < SMALL_COHORT:
    L += ['', f'> **MECHANICS TEST ONLY.** The smallest group has {min(inp0["groups"].values())} sample(s) (fewer than {SMALL_COHORT}). With so few samples the DESeq2 results, the candidate-rule results (including the leave-one-out robustness rule), the pathway analysis and the LASSO model in this report are a mechanics check only, NOT a biomarker result. Do not interpret p-values, gene lists, pathways, the model or its performance.']
L += ['', '## 1. Input', '']
if inp:
    L += [f'{inp["n_samples"]} samples ({", ".join(f"{k} {v}" for k, v in inp["groups"].items())}); data: {inp["data"]}; '
          f'scanned folder entries: {inp["n_scanned_entries"]}, problem files: {inp["n_problem_files"]}.', '',
          '**Group label mapping used** (conf/input.config `group_aliases`; case-insensitive, whole string):', '',
          md_table(rows(J('input', 'group_mapping.tsv')), ['label_in_metadata', 'canonical_group', 'n_samples'])]
    if inp['unknown_samples']:
        L += [f'Samples scored without a label (group = unknown): {", ".join(inp["unknown_samples"])}', '']
    for n in inp['notes']:
        L += [f'- {n}']
    L += ['batch and study: from `batch_or_study`, unless separate `batch` / `study` columns exist (owner decision). '
          'patient_id, where given, is carried and not used in any score.', '']
else:
    L += ['The run stopped at the input check; see `input_scan.log` and `input/problem_files.tsv`.', '']
if os.path.exists(J('input', 'qc_excluded.tsv')):
    qx = rows(J('input', 'qc_excluded.tsv'))
    L += [f'Samples excluded by the sample QC (status FAIL): {len(qx)}' + ('' if not qx else ' - ' + '; '.join(f"{r['sample']} ({r['group']}): {r['reason']}" for r in qx)), '']
if os.path.exists(J('input', 'excluded_after_quant.tsv')):
    ex = rows(J('input', 'excluded_after_quant.tsv'))
    L += [f'Samples excluded after QUANT (QC FAIL): {len(ex)}', '']

ST = []
if info.get('stop_after'):
    ST += [f"- Run stopped after stage `{info['stop_after']}` (--stop_after): later stages were not run. Continue with the same command without --stop_after and with -resume."]
for st in ('quant', 'discover'):
    e = stage_exit(st)
    if e is not None:
        ST += [f'- Nextflow stage `{st}`: exit {e} (logs in `nf_{st}/`)']
_tr = sorted(glob.glob(J('logs', 'trace_*.tsv')) + glob.glob(J('trace*.tsv')))
if _tr:                                    # one-workflow route: Nextflow trace of the run (if it is part of the run folder)
    _tn = {}
    for _r in rows(_tr[-1]):
        _k = (_r.get('name') or '').split(' (')[0]
        _tn.setdefault(_k, []).append(_r.get('status', ''))
    ST += [f'- Nextflow trace `{os.path.relpath(_tr[-1], R)}`: ' + '; '.join(f'{k} x{len(v)} ({", ".join(sorted(set(v)))})' for k, v in _tn.items() if k)]
elif not ST:                               # no trace and no nf_<stage> folder: list the stage output folders that exist (this report was built from them)
    _sf = [(d, len([f for f in glob.glob(J('pipeline', d, '**', '*'), recursive=True) if os.path.isfile(f)])) for d in ('quant', 'analysis', 'discover', 'apply') if os.path.isdir(J('pipeline', d))]
    if _sf:
        ST += ['- Stage outputs present in this run (no Nextflow trace was available to the report step): ' + '; '.join(f'`pipeline/{d}/` ({n} files)' for d, n in _sf)]
if os.path.exists(J('pipeline', 'apply', 'apply.log')):
    ST += [f'- locked-model scoring (`scripts/apply_locked_model.R`, lock/diagnostic_v1, hash-checked): see `pipeline/apply/apply.log`']
if ST:
    L += ['## 2. Stages', ''] + ST + ['']

if MODE == 'discover' and os.path.exists(J('pipeline', 'discover', '03_model', 'model_spec.txt')):
    spec = open(J('pipeline', 'discover', '03_model', 'model_spec.txt')).read().splitlines()
    cut = next((s for s in spec if s.startswith('cutoff = ')), 'cutoff: not found')
    mg = rows(J('pipeline', 'discover', '03_model', 'model_genes.tsv'))
    rr_path = J('pipeline', 'discover', '02_candidates', 'robustness_refits.tsv')
    rob = []
    _ds = J('pipeline', 'discover', '01_de', 'deseq2_summary.json')       # the DESeq2 design actually used (batch term only when there is more than one batch)
    if os.path.exists(_ds):
        _dz = (json.load(open(_ds)) or {}).get('design')
        if _dz:
            rob += [f'**DESeq2 design:** `{_dz}`' + (' (one batch only: no batch term is possible)' if _dz.strip() == '~ group' else ' (batch is adjusted for)'), '']
    if os.path.exists(rr_path):
        rr = rows(rr_path)
        if rr:
            ran = [r for r in rr if r['status'] == 'run']; skp = [r for r in rr if r['status'] != 'run']
            rob = [f'**Robustness refits (rule 3b):** {len(ran)} of {len(rr)} leave-one-out refits ran, {len(skp)} skipped.']
            if skp:
                rob += ['Skipped refits (the design `batch + group` needs every batch to keep both groups): ' + '; '.join(f"dropping {r['refit']}: {r['note'].replace('skipped: ', '')}" for r in skp) + '. The rule 3b threshold is 90 % of the refits that ran (rounded up); with fewer than 10 refits run, rule 3b is NOT EVALUABLE and no candidate is declared.']
            rob += ['']
    _ss, _nt = J('pipeline', 'discover', '00_input', 'sample_sheet.discover.tsv'), J('pipeline', 'discover', '02_candidates', 'candidate_rules_method_note.txt')
    if os.path.exists(_ss) and os.path.exists(_nt) and 'rule 3b: loo' in open(_nt).read() and 'leave-one-PATIENT-out' not in open(_nt).read():
        _pid = [r.get('patient_id') for r in rows(_ss) if r.get('patient_id')]
        if len(set(_pid)) < len(_pid):
            rob += ['Repeated patient_id found; folds and refits are sample-level (the locked method). Use --patient_level_cv for patient-level.', '']
    notfit = os.path.exists(J('pipeline', 'discover', '03_model', 'model_spec.txt')) and 'NOT FITTED' in open(J('pipeline', 'discover', '03_model', 'model_spec.txt')).read()
    _why = ''
    if os.path.exists(_nt):               # the real reason, as written by the candidate-rules step (R5)
        _why = next((l.strip() for l in open(_nt).read().splitlines() if 'RULE 3B NOT EVALUABLE' in l), '')
    L += ['## 3. Discovery result (development fit)', '', '**Candidate funnel:**', '', md_table(rows(J('pipeline', 'discover', '02_candidates', 'funnel_counts.tsv'))), *rob,
          (('**NO MODEL WAS FITTED:** ' + (f'{_why}. ' if _why else 'fewer than 2 candidate genes passed all applied rules (see the funnel above), so the LASSO did not run. ') + 'This is not a negative biomarker result.') if notfit else f'**Candidate model** `model/diagnostic_v2_candidate_{RID}/` (NOT locked; not in lock/): {len(mg)} genes, {cut}'), '',
          *([md_table(mg, [c for c in ('gene_id', 'gene_name', 'sign', 'lasso_coef_lambda1se') if mg and c in mg[0]])] if not notfit else []),
          '**Development performance (OPTIMISTIC, not validation):**', '', md_table(rows(J('pipeline', 'discover', '03_model', 'development_performance.tsv'))),
          *(['**Sensitivity and specificity (Clopper-Pearson 95% CI; development, OPTIMISTIC):**', '', md_table(rows(J('pipeline', 'discover', '03_model', 'development_sens_spec.tsv')))] if os.path.exists(J('pipeline', 'discover', '03_model', 'development_sens_spec.tsv')) else []),
          *(['Pathway analysis (over-representation and gene-set enrichment, GO BP / Hallmark) ran on the differential-expression results of all QC-pass samples; see the Plots section and `pipeline/analysis/07_pathway/`. It is exploratory and is not used by the candidate rules or the model.', ''] if any(st == 'analysis/07_pathway' for st in plot_src.values()) else [])]
if MODE == 'apply' and os.path.exists(J('pipeline', 'apply', 'scores', 'scores.tsv')):
    L += ['## 3. Locked-model scores (lock/diagnostic_v1: genes, signs and cutoff from the lock; the lock was built on colon tumour vs normal)', '', md_table(rows(J('pipeline', 'apply', 'scores', 'scores.tsv')))]
    if os.path.exists(J('pipeline', 'apply', 'scores', 'evaluation.tsv')):
        L += ['**Evaluation against the metadata labels (report-only):**', '', md_table(rows(J('pipeline', 'apply', 'scores', 'evaluation.tsv')))]
        if os.path.exists(J('pipeline', 'apply', 'scores', 'evaluation_sens_spec.tsv')):
            L += ['**Sensitivity and specificity (Clopper-Pearson 95% CI):**', '', md_table(rows(J('pipeline', 'apply', 'scores', 'evaluation_sens_spec.tsv')))]
    L += ['Cohort outlier flag: not computed for small cohorts (single-sample baseline DEFERRED).', '']
PLOT_GROUPS = [('quant/05_cohort_qc', 'Sample QC (all quantified samples)'), ('analysis/06_dge', 'Differential expression, all QC-pass samples'),
               ('analysis/07_pathway', 'Pathway analysis (over-representation and gene-set enrichment)'), ('discover/01_de', 'Differential expression, discovery training samples')]
CAPTION = [('volcano', 'Volcano plot'), ('_MA', 'MA plot'), ('top50_heatmap', 'Heatmap of the top 50 genes'), ('pca_group_study', 'PCA of the samples (colour = group, shape = study)'),
           ('sample_correlation_heatmap', 'Sample correlation heatmap'), ('library_size_barplot', 'Library sizes'), ('dotplot', 'Dot plot of enriched pathways'), ('_top', 'Top enriched gene sets')]
cap = lambda n: next((c for k, c in CAPTION if k in n), n)
pl = ['', '']
had = False
# R6: the discovery stage re-runs DESeq2 on its training samples. When these are all the QC-pass samples (no hold-out), say so instead of a second look-alike section.
_dz2 = J('pipeline', 'discover', '01_de', 'deseq2_summary.json')
_nall = (inp['n_samples'] - len(rows(J('input', 'qc_excluded.tsv'))) - len(rows(J('input', 'excluded_after_quant.tsv')))
         if inp and os.path.exists(J('input', 'qc_excluded.tsv')) and os.path.exists(J('input', 'excluded_after_quant.tsv')) else inp['n_samples'] if inp else None)
SAME_DE = bool(_nall and os.path.exists(_dz2) and (json.load(open(_dz2)) or {}).get('n_samples') == _nall
               and any(st == 'analysis/06_dge' for st in plot_src.values()) and any(st == 'discover/01_de' for st in plot_src.values()))
for stage, title in PLOT_GROUPS:
    names = sorted(n for n, st in plot_src.items() if st == stage)
    if SAME_DE and stage == 'discover/01_de':
        dup = [n for n in names if n.startswith('01_de_') and plot_src.get(n[len('01_de_'):]) == 'analysis/06_dge']
        names = [n for n in names if n not in dup]
        title = 'Differential expression, discovery training samples'
        pl += ['The discovery stage re-ran differential expression on its training samples, which are the same samples as above (no hold-out). Its plots for the same route are the files `plots/01_de_*` and are not repeated here' + ('; plots that exist only in the discovery stage follow.' if names else '.'), '']
        had = True
        if not names:
            continue
        pl += ['**Discovery stage only (other route)**', ''] + [f'![{cap(n)}](../plots/{n})' for n in names] + ['']
        continue
    if names:
        had = True
        pl += [f'**{title}**', ''] + [f'![{cap(n)}](../plots/{n})' for n in names] + ['']
if had:
    L += ['## Plots (files are also in `plots/`)', ''] + pl[2:]
# ---- optional: internal validation written by the discovery stage (only if the folder exists)
_vd = J('pipeline', 'discover', '04_validation')
if os.path.isdir(_vd):
    L += ['## Internal validation (pipeline/discover/04_validation/)', '']
    if os.path.exists(os.path.join(_vd, 'internal_validation_summary.tsv')):
        L += ['**Internal validation summary** (same samples as the development fit; not independent validation):', '',
              md_table(rows(os.path.join(_vd, 'internal_validation_summary.tsv')), ['metric', 'value', 'ci_low', 'ci_high', 'label']), '']
    if os.path.exists(os.path.join(_vd, 'selection_frequency.tsv')):
        L += ['**Selection frequency of genes** (how often each gene was selected on resampled data):', '',
              md_table(rows(os.path.join(_vd, 'selection_frequency.tsv')), ['gene_id', 'symbol', 'direction', 'frequency', 'in_final_model']), '']
    if os.path.exists(os.path.join(_vd, 'INTERNAL_VALIDATION.md')):
        L += ['**INTERNAL_VALIDATION.md:**', ''] + [f'**{l.lstrip("# ")}**' if l.startswith('#') else l for l in open(os.path.join(_vd, 'INTERNAL_VALIDATION.md')).read().splitlines()] + ['']
_ns = os.path.join(_vd, 'nested_validation_summary.tsv')
if os.path.exists(_ns):
    L += ['## Nested internal validation (candidate screen repeated inside every training fold)', '',
          'Rules 1-2, LASSO, score and cutoff are all re-done on the training samples of each fold and the held-out fold is scored with that model, so the gene selection cannot see the samples it is scored on. Rules 3a and 3b are not run, so this estimates the simplified procedure, not the full locked one. It is internal, not independent validation. The permutation null shows what AUC the same procedure reaches on shuffled labels.', '',
          md_table(rows(_ns), ['metric', 'value', 'ci_low', 'ci_high', 'label']), '']
# ---- warnings written by input scan / cohort QC (lines that begin with WARNING); input notes are already listed in section 1
_seen = set(inp['notes']) if inp else set()
_w = []
for _f in sorted(glob.glob(J('pipeline', 'quant', '05_cohort_qc', '*'))):
    if os.path.isfile(_f) and _f.endswith(('.txt', '.md', '.tsv', '.log')):
        _w += [l.strip() for l in open(_f, errors='replace').read().splitlines() if l.startswith('WARNING')]
    elif os.path.isfile(_f) and _f.endswith('.json'):
        try:
            _j = json.load(open(_f))
        except ValueError:
            continue
        _w += [n for k in ('notes', 'warnings') for n in (_j.get(k) or []) if isinstance(_j, dict) and isinstance(n, str) and n.startswith('WARNING')]
_w = [w for i, w in enumerate(_w) if w not in _seen and w not in _w[:i]]
if _w:
    L += ['## Warnings from the cohort QC', ''] + [f'- {w}' for w in _w] + ['']
_rqp = J('pipeline', 'quant', '05_cohort_qc', 'rna_qc_metrics.tsv')
_rq = rows(_rqp) if os.path.exists(_rqp) else []
if _rq:
    _fl = [r for r in _rq if r.get('flags') and not r['flags'].startswith('not enough')]
    L += ['## RNA quality metrics (flag only)', '',
          'Per-sample metrics from RSeQC and the featureCounts file: exonic and intronic read share, position-based duplication, transcript integrity number (TIN, median over a seeded random subset of transcripts), 3\'/5\' gene-body coverage ratio (our own definition, not Picard\'s), rRNA and mitochondrial share of assigned reads (rRNA is a lower bound: only annotated rRNA genes are counted, rDNA arrays and multi-mapped reads are not). '
          'No published pass/fail cut-offs exist for these metrics, so a sample is only flagged when it is an outlier relative to the rest of this cohort (robust z-score beyond the cohort-outlier cut-off); with fewer than 5 samples nothing is flagged. Flags never exclude a sample or change a QC status. When a cohort mixes studies or library preparations the metrics differ systematically between them (in the 6-sample demo the two samples of one study show strong 3\' bias and a high intronic share), so flags then mostly show study differences: compare within a study.', '',
          (f'**{len(_fl)} of {len(_rq)} sample(s) flagged:** ' + '; '.join(f"{r['srr']}: {r['flags']}" for r in _fl)) if _fl else ('No sample is an outlier on these metrics.' if len(_rq) >= 5 else f'Only {len(_rq)} samples: metrics are listed, no cohort flags are computed.'), '',
          md_table(_rq, [c for c in ('srr', 'exonic_pct', 'intronic_pct', 'dup_rate_pct', 'tin_median', 'bias_3p_5p', 'rrna_pct', 'mito_pct') if c in _rq[0]]), '']
_scp = J('pipeline', 'quant', '05_cohort_qc', 'strandedness_check.tsv')
_sc = rows(_scp) if os.path.exists(_scp) else []
if _sc:
    _bad = [r for r in _sc if r.get('agree') == 'NO']
    _und = [r for r in _sc if r.get('agree') == 'undetermined']
    L += ['## Strandedness cross-check (flag only)', '',
          'featureCounts used the strandedness read from the Salmon library type. RSeQC `infer_experiment.py` gives a second, independent estimate from the alignment. Nothing below changes a count, a QC status or an exclusion.', '',
          (f'**WARNING: {len(_bad)} of {len(_sc)} sample(s) disagree (' + ', '.join(r['srr'] for r in _bad) + '). Check the library type: gene counts of these samples may be wrong.**') if _bad else f'All {len(_sc) - len(_und)} determinable sample(s) agree with the Salmon call.',
          *([f'{len(_und)} sample(s) undetermined by RSeQC (ambiguous fractions; can also indicate genomic DNA contamination): ' + ', '.join(r['srr'] for r in _und)] if _und else []), '',
          md_table(_sc, [c for c in ('srr', 'rseqc_forward', 'rseqc_reverse', 'rseqc_failed', 'rseqc_call', 'salmon_format', 'featurecounts_s', 'agree') if c in _sc[0]]), '']
# ---- methods and versions (from files that already exist)
M = ['## Methods and versions', '']
_tv = {}
_gr = (json.load(open(_mq0)) if os.path.exists(_mq0 := J('pipeline', 'quant', 'manifest_quant.json')) else {})
_tv = _gr.get('tool_versions') or {}
if _tv:
    M += ['Tool versions as recorded by the quantification stage (`manifest_quant.json`):', '', '| tool | version |', '|--|--|'] + [f'| {k} | {v} |' for k, v in _tv.items()] + ['']
else:
    M += ['Tool versions: not recorded in this run (no `tool_versions` in `manifest_quant.json`); see `tool_versions.tsv`.', '']
_gv = (_gr.get('parameters') or {}).get('gencode_release')
M += [f'**Annotation:** GENCODE release {_gv}.' if _gv else '**Annotation:** GENCODE release not recorded in this run.']
_rf = _gr.get('references') or {}
if _rf:
    def _ref(k, v):
        if isinstance(v, dict) and 'path' in v:
            return f'- {k}: `{os.path.basename(str(v["path"]))}`, {v.get("size")} bytes' + (f', md5 `{v["md5"]}`' if v.get('md5') else '')
        if isinstance(v, dict):                 # an index: several files, md5 only where the quantification recorded one
            return f'- {k}: ' + '; '.join(f'`{n}` ({f.get("size")} bytes' + (f', md5 `{f["md5"]}`' if f.get('md5') else '') + ')' for n, f in v.items() if isinstance(f, dict))
        return f'- {k}: {v}'
    M += ['', 'Reference files as recorded by the quantification (md5 where recorded; large index files are listed by size only):', ''] + [_ref(k, v) for k, v in _rf.items()]
_th = os.path.join(P, 'conf', 'thresholds.config')
M += ['', '**Thresholds file** `conf/thresholds.config` sha256: ' + (f'`{sha(_th)}`' if os.path.exists(_th) else 'not found') + ' (values marked PROVISIONAL in it have no external source).', '']
L += M
# ---- limitations (fixed wording; the small-cohort line only when it applies)
L += ['## Limitations and recommended next validation steps', '',
      '- Research use only. This output is not a clinical test and is not clinically validated.',
      '- Performance figures from the samples used to choose the genes and cutoff are in-sample and optimistic.',
      '- Candidate genes are hypotheses. They need independent validation in samples that were not used to build them.']
if not info.get('counts'):
    L += ['- The Salmon route is not exactly repeatable: repeated runs on the same reads can give slightly different Salmon counts and DESeq2 statistics, so borderline candidates can change. The featureCounts route is repeatable (identical counts in the one repeat test of this pipeline).']
if inp and inp.get('groups'):
    _n, _g = inp['n_samples'], min(inp['groups'].values())
    if _n < 10 or _g < 6:
        L += [f'- Small cohort: {_n} samples, smallest group {_g}. Fewer than 10 samples or fewer than 6 per group gives low power and unstable gene lists; treat all gene-level results as exploratory.']
L += ['- Next step: if a model was fitted, freeze its genes, signs and cutoff before looking at new data, then apply them unchanged to an independent cohort (different study and site, none of whose samples were used here) and report sensitivity, specificity and their confidence intervals. The locked model of this repository is specific to colon tumour versus normal and is not a result of this run.', '']
# ---- reporting statement (TRIPOD+AI / REMARK / STARD-oriented; fixed wording, says what this run does NOT establish)
L += ['## Reporting statement', '',
      '- **Intended use:** research use only; bulk RNA-seq of human tissue, case versus reference groups as given in the metadata. Sample type (fresh-frozen or FFPE), tissue, site and stage are not checked by the pipeline. The locked 9-gene model of this repository (apply mode) is specific to colon tumour versus normal.',
      '- **Reference standard:** the group labels are taken from `metadata.csv` as given. The pipeline does not verify them against histology; the only identity check is the flag-only sex check in the cohort QC.',
      '- **Blinding and partitioning:** the pipeline does not know which samples are independent. Every development number is computed on the samples used to choose genes and cutoff; only an apply run on a cohort that was never used (`--cohort validation`) with a model frozen beforehand is an evaluation.',
      '- **Candidates considered and selection steps:** the candidate funnel above lists every gene entering each rule; the rules, thresholds and seeds are in `params.json` and `conf/thresholds.config`. Candidate screening is not repeated inside the internal cross-validation.',
      '- **Model specification:** genes, signs, cutoff and training-sample list are in `pipeline/discover/03_model/` (or `lock/diagnostic_v1/` for a locked run); `manifest.json` holds the sha256 of every file, and `run_info.json` the code commit and the uncommitted files at the start.',
      '- **Not assessed here:** calibration of the score, missing-data handling, performance by subgroup (site, sex, stage), fairness and applicability to other populations. Confidence intervals above are exact (Clopper-Pearson) and are wide for small groups.', '']
L += ['## 4. Files', '', '`metadata.used.csv` (copy of the metadata used), `input/` (scan listing, resolved metadata, sheets), `pipeline/` (all outputs), '
      f'`tables/`, `plots/`, `tool_versions.tsv`, `manifest.json` (sha256 of inputs, config, lock and report files). Links made now: {", ".join(made) or "none"}.', '']
text = '\n'.join(L) + '\n'


def to_html(md):
    out, in_tab, in_list = [], False, False
    for line in md.splitlines():
        esc = html.escape(line)
        esc = __import__('re').sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', esc)
        esc = __import__('re').sub(r'`(.+?)`', r'<code>\1</code>', esc)
        mimg = __import__('re').fullmatch(r'!\[(.*?)\]\((.*?)\)', line.strip())
        if mimg:
            if in_tab:
                out.append('</table>'); in_tab = False
            out.append(f'<figure style="display:inline-block;margin:.3em;max-width:48%"><img src="{html.escape(mimg.group(2))}" alt="{html.escape(mimg.group(1))}" style="max-width:100%"><figcaption>{html.escape(mimg.group(1))}</figcaption></figure>')
            continue
        if line.startswith('|'):
            cells = [c.strip() for c in esc.strip('|').split('|')]
            if not in_tab:
                out.append('<table><tr>' + ''.join(f'<th>{c}</th>' for c in cells) + '</tr>'); in_tab = True
            elif not set(line) <= set('|-'):
                out.append('<tr>' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>')
            continue
        if in_tab:
            out.append('</table>'); in_tab = False
        if line.startswith('- '):
            if not in_list:
                out.append('<ul>'); in_list = True
            out.append(f'<li>{esc[2:]}</li>'); continue
        if in_list:
            out.append('</ul>'); in_list = False
        if line.startswith('#'):
            n = len(line) - len(line.lstrip('#'))
            out.append(f'<h{n}>{esc[n:].strip()}</h{n}>')
        elif line.startswith('&gt; '):
            out.append(f'<blockquote>{esc[5:]}</blockquote>')
        elif line.strip():
            out.append(f'<p>{esc}</p>')
    if in_tab:
        out.append('</table>')
    css = ('body{font-family:sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;line-height:1.4}table{border-collapse:collapse;margin:.5em 0}'
           'td,th{border:1px solid #bbb;padding:2px 6px;font-size:90%}blockquote{background:#fff3cd;padding:.5em;margin:0}code{background:#eee}')
    return f'<!doctype html><meta charset="utf-8"><title>{html.escape(RID)}</title><style>{css}</style>\n' + '\n'.join(out) + '\n'


done = att[-1].get('status') == 'completed'
rep = J('report') if done and not os.path.exists(J('report')) else J(f'report_{TS}' if done else f'report_partial_{TS}')
os.makedirs(rep)
open(os.path.join(rep, 'report.md'), 'w').write(text)
open(os.path.join(rep, 'report.html'), 'w').write(to_html(text))

# ---- 3. manifest
files = {}
for d in ('input', 'tables', 'model', os.path.basename(rep)):
    for p in sorted(glob.glob(J(d, '**', '*'), recursive=True)):
        if os.path.isfile(p) and not os.path.islink(p):
            files[os.path.relpath(p, R)] = sha(p)
def nf_version():
    """Nextflow version from env/bin/nextflow of this checkout; 'unknown' when that binary is not there (e.g. a worktree without env/)."""
    try:
        return subprocess.run([os.path.join(P, 'env/bin/nextflow'), '-version'], capture_output=True, text=True).stdout.split('version')[-1].split('\n')[0].strip()
    except OSError:
        return 'unknown'


ins = {k: v for k, v in (('metadata.used.csv', J('metadata.used.csv')), ('counts', info.get('counts'))) if v and os.path.exists(v)}
man = dict(run_id=RID, mode=MODE, written=TS, git_commit=info['git_commit'], uncommitted_files=info['uncommitted_files'],
           command=info['command'], inputs={k: dict(path=v, sha256=sha(v)) for k, v in ins.items()},
           config={c: sha(os.path.join(P, c)) for c in ('conf/input.config', 'conf/thresholds.config', 'config/candidate_rules.yaml', 'config/questions.yaml')},
           lock={l.split()[1]: l.split()[0] for l in open(os.path.join(P, 'lock/diagnostic_v1/sha256sums_final.txt'))},
           nextflow_version=nf_version(),
           files=files)
_mq = J('pipeline', 'quant', 'manifest_quant.json')       # reference files as recorded by the quantification (md5 and sizes; no new hashing)
if os.path.exists(_mq):
    _m = json.load(open(_mq))
    man['reference'] = dict(gencode_release=(_m.get('parameters') or {}).get('gencode_release'), files=_m.get('references'),
                            note='copied from manifest_quant.json: md5 and sizes (md5 only; the whole STAR/Salmon index is not hashed)')
mp = J('manifest.json') if not os.path.exists(J('manifest.json')) else J(f'manifest_{TS}.json')
json.dump(man, open(mp, 'w'), indent=2)
print(f'report: {os.path.join(rep, "report.md")} (+ .html); manifest: {mp}')
