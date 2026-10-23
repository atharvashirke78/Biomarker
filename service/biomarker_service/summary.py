"""Read-only summaries of a job's output folder for the UI: sample QC counts, candidate funnel, model, development performance, figures.
Everything shown is copied from the pipeline's own files; the limitation texts are fixed and never say that anything is validated."""
import csv
import os

from . import db

CAPTIONS = [('volcano', 'Volcano plot'), ('_MA', 'MA plot'), ('top50_heatmap', 'Heatmap of the top 50 genes'),
            ('pca_group_study', 'PCA of the samples (colour = group, shape = study)'), ('sample_correlation_heatmap', 'Sample correlation heatmap'),
            ('library_size_barplot', 'Library sizes'), ('dotplot', 'Dot plot of enriched pathways'), ('_top', 'Top enriched gene sets')]
LIMITATIONS = [
    'Automated summary of the pipeline outputs. Interpretation is a separate step done by the research team.',
    'Development performance is computed on the same samples that were used to choose the genes and the cutoff: it is optimistic and is not validation.',
    'Candidate biomarkers are hypotheses. They need independent validation on samples that were not used here before any diagnostic claim.',
    'Research use only. Not a clinical or diagnostic device.']


def out_dir(job_id):
    return db.paths(job_id)['output']


def read_tsv(job_id, rel):
    path = os.path.join(out_dir(job_id), *rel.split('/'))
    try:
        with open(path, newline='', encoding='utf-8', errors='replace') as f:
            return list(csv.DictReader(f, delimiter='\t'))
    except OSError:
        return None


def read_text(job_id, rel):
    try:
        with open(os.path.join(out_dir(job_id), *rel.split('/')), errors='replace') as f:
            return f.read()
    except OSError:
        return None


def figures(job_id):
    d = os.path.join(out_dir(job_id), '10_report', 'plots')
    if not os.path.isdir(d):
        return []
    res = []
    for n in sorted(os.listdir(d)):
        if not n.lower().endswith('.png') or n.startswith('.'):
            continue
        cap = next((c for k, c in CAPTIONS if k in n), n)
        if n.startswith('01_de_'):
            cap += ' (discovery training samples)'
        res.append(dict(name=n, caption=cap, path='10_report/plots/' + n))
    return res


def summary(job_id):
    s = dict(limitations=LIMITATIONS)
    qc = read_tsv(job_id, '04_counts/sample_qc.tsv')
    if qc is not None:
        cnt = {}
        for r in qc:
            cnt[r.get('status', '?')] = cnt.get(r.get('status', '?'), 0) + 1
        s['qc'] = dict(n_samples=len(qc), counts=cnt,
                       flagged=[dict(sample=r['srr'], status=r['status'], reason=r.get('reason', ''), group=r.get('group', '')) for r in qc if r.get('status') != 'PASS'])
    fun = read_tsv(job_id, '09_discover/02_candidates/funnel_counts.tsv')
    if fun is not None:
        s['funnel'] = [dict(step=r['step'], all_genes=r['all_genes'], protein_coding=r.get('protein_coding', '')) for r in fun]
    cand = read_tsv(job_id, '09_discover/02_candidates/candidate_table.tsv')
    if cand is not None:
        passed = [r for r in cand if r.get('all_rules_pass') == 'TRUE']
        s['candidates'] = dict(n_table=len(cand), n_pass_all_rules=len(passed),
                               top=[dict(symbol=r.get('symbol', ''), gene_id=r['gene_id'], log2fc_shrunk=r.get('log2FC_shrunk_fc', ''), padj=r.get('padj_fc', ''),
                                         loo_passes=r.get('loo_n_pass', '')) for r in passed[:15]])
    spec = read_text(job_id, '09_discover/03_model/model_spec.txt')
    genes = read_tsv(job_id, '09_discover/03_model/model_genes.tsv')
    if spec is not None:
        cut = next((l.split('=', 1)[1].split('(')[0].strip() for l in spec.splitlines() if l.startswith('cutoff = ')), None)
        fitted = 'NOT FITTED' not in spec and bool(genes)
        s['model'] = dict(fitted=fitted, cutoff=cut, genes=[dict(symbol=g.get('symbol', ''), gene_id=g['gene_id'], direction=g.get('set', ''), coef=g.get('lasso_coef_lambda1se', '')) for g in (genes or [])],
                          note=None if fitted else next((l for l in spec.splitlines() if l.startswith('NO MODEL FITTED')), 'no model was fitted'))
    perf = read_tsv(job_id, '09_discover/03_model/development_performance.tsv')
    sens = read_tsv(job_id, '09_discover/03_model/development_sens_spec.tsv')
    if perf and 'AUC' in perf[0]:
        s['performance'] = dict(rows=perf, sens_spec=sens if sens and 'sensitivity' in sens[0] else None)
    return s
