#!/usr/bin/env python3
"""Regression checks for the 19-sample DEVELOPMENT outputs (read-only).

Compares result files against the documented expected values in expected_19_development.json.
Prints only counts and PASS/FAIL, never expression values. Development data only: nothing here is validation.

usage (from the project root, pyqc env):
  env/prefix/pyqc/bin/python3 -I tests/regression/check_19_development.py --out <summary.tsv>
      [--de-dir results_rc1prep_19] [--candidate-table results_rc1prep_19/rule3b_rerun18_cfg/candidate_table.tsv]
      [--final-candidates results_rc1prep_19/final_candidates/final_candidates_18of19.tsv]
      [--reactome-dir results_rc1prep_19/reactome] [--scores <locked-model scores.tsv from apply_locked_model.R>]
      [--rule4 <rule4 per-gene table from rule4_reproduction.R, self-consistency run on the 19-sample fit>]
Exit code 0 if every check passes, 1 otherwise.
"""
import argparse
import hashlib
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []


def check(name, observed, expected, ok=None, source=''):
    if ok is None:
        ok = observed == expected
    rows.append(dict(check=name, observed=observed, expected=expected, status='PASS' if ok else 'FAIL', source=source))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def de_counts(path):
    d = pd.read_csv(path, sep='\t')
    r1 = d.padj.lt(0.05) & d.log2FoldChange_shrunk.abs().ge(1)
    return dict(genes_after_prefilter=len(d), tested_nonNA_pvalue=int(d.pvalue.notna().sum()),
                **{'padj_lt_0.05': int(d.padj.lt(0.05).sum())},
                rule1_total=int(r1.sum()), rule1_up=int((r1 & (d.log2FoldChange_shrunk > 0)).sum()),
                rule1_down=int((r1 & (d.log2FoldChange_shrunk < 0)).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--expected', default=os.path.join(HERE, 'expected_19_development.json'))
    ap.add_argument('--de-dir', default='results_rc1prep_19')
    ap.add_argument('--candidate-table', default='results_rc1prep_19/rule3b_rerun18_cfg/candidate_table.tsv')
    ap.add_argument('--final-candidates', default='results_rc1prep_19/final_candidates/final_candidates_18of19.tsv')
    ap.add_argument('--reactome-dir', default='results_rc1prep_19/reactome')
    ap.add_argument('--sample-sheet', default='results_rc1prep_19/sample_sheet_rc1prep.tsv')
    ap.add_argument('--patient-map', default='results/09_verification/sample_map_with_patient.tsv')
    ap.add_argument('--scores', default=None)
    ap.add_argument('--rule4', default=None)
    ap.add_argument('--featurecounts-only', action='store_true', help='check only the featureCounts DE route (the -entry ANALYSIS output has no Salmon route)')
    a = ap.parse_args()
    E = json.load(open(a.expected))

    # 0. lock integrity (sha256sums_final.txt paths are relative to the project root)
    for line in open('lock/diagnostic_v1/sha256sums_final.txt'):
        h, p = line.split()
        check(f'lock sha256 {p}', 'match' if sha256(p) == h else 'MISMATCH', 'match', source='lock/diagnostic_v1/sha256sums_final.txt')

    # 1. cohort
    ss = pd.read_csv(a.sample_sheet, sep='\t')
    c = E['cohort']
    check('cohort n_samples', len(ss), c['n_samples'], source=a.sample_sheet)
    check('cohort n_normal', int((ss.group == 'normal').sum()), c['n_normal'], source=a.sample_sheet)
    check('cohort n_tumour', int((ss.group == 'tumour').sum()), c['n_tumour'], source=a.sample_sheet)
    for st, g in c['per_study'].items():
        for grp, n in g.items():
            check(f'cohort {st} {grp}', int(((ss.study == st) & (ss.group == grp)).sum()), n, source=a.sample_sheet)
    pm = pd.read_csv(a.patient_map, sep='\t')
    pm = pm[pm.srr.isin(ss.sample_id)]
    check('cohort samples in patient map', len(pm), c['n_samples'], source=a.patient_map)
    check('cohort n_patients', pm.patient_id.nunique(), c['n_patients'], source=a.patient_map)
    pairs = pm.groupby('patient_id').group.nunique().eq(2).sum()
    check('cohort tumour/normal pairs', int(pairs), c['n_tumour_normal_pairs'], source=a.patient_map)

    # 2. DE per route
    routes = [('de_featurecounts', 'deseq2_results_all.tsv')] + ([] if a.featurecounts_only else [('de_salmon', 'deseq2_salmon_results_all.tsv')])
    for key, fn in routes:
        p = os.path.join(a.de_dir, fn)
        if not os.path.exists(p):
            check(f'{key} file present', 'missing', 'present', source=p)
            continue
        obs = de_counts(p)
        for k, v in E[key].items():
            if not k.startswith('_'):
                check(f'{key} {k}', obs[k], v, source=p)

    # 3. funnel and final candidates
    if os.path.exists(a.candidate_table):
        ct = pd.read_csv(a.candidate_table, sep='\t')
        t = lambda col: ct[col].astype(str).eq('True')
        f = E['funnel']
        r1, r2, r3a, fin = t('rule1_pass'), t('rule2_pass'), t('rule3a_pass'), t('all_rules_pass')
        check('funnel rule1 either route', int(r1.sum()), f['rule1_either'], source=a.candidate_table)
        check('funnel + rule2', int((r1 & r2).sum()), f['rule2'], source=a.candidate_table)
        check('funnel + rule3a', int((r1 & r2 & r3a).sum()), f['rule3a'], source=a.candidate_table)
        check('funnel + rule3b (final)', int(fin.sum()), f['rule3b_final'], source=a.candidate_table)
        check('funnel final protein_coding', int((fin & ct.gene_biotype.eq('protein_coding')).sum()), f['rule3b_final_protein_coding'], source=a.candidate_table)
        check('funnel baseMean>=50 tier', int((fin & t('baseMean_ge_50_both')).sum()), f['basemean50_tier'], source=a.candidate_table)
        fc = pd.read_csv(a.final_candidates, sep='\t')
        check('final list == all_rules_pass genes', set(fc.gene_id) == set(ct.gene_id[fin]), True, source=a.final_candidates)
    else:
        check('candidate table present', 'missing', 'present', source=a.candidate_table)
    fc = pd.read_csv(a.final_candidates, sep='\t')
    e = E['final_candidates']
    check('final candidates n', len(fc), e['n'], source=a.final_candidates)
    check('final candidates up', int(fc.direction_fc.eq('up').sum()), e['up'], source=a.final_candidates)
    check('final candidates down', int(fc.direction_fc.eq('down').sum()), e['down'], source=a.final_candidates)
    check('final candidates protein_coding', int(fc.gene_biotype.eq('protein_coding').sum()), e['protein_coding'], source=a.final_candidates)

    # 4. Reactome
    rs = pd.read_csv(os.path.join(a.reactome_dir, 'analysis_summary.tsv'), sep='\t').set_index('analysis')
    r = E['reactome']
    for k, lab in [('ora_up_significant', 'ORA up'), ('ora_down_significant', 'ORA down'), ('ora_combined_significant', 'ORA combined'), ('gsea_significant', 'GSEA')]:
        check(f'reactome {lab} significant', int(rs.loc[lab, 'sets_significant']), r[k], source=a.reactome_dir)

    # 5. locked model artefacts (and, if given, scores from apply_locked_model.R)
    m = E['locked_model']
    g = pd.read_csv('lock/diagnostic_v1/final_model_genes.tsv', sep='\t')
    check('locked model n_genes', len(g), m['n_genes'], source='lock/diagnostic_v1/final_model_genes.tsv')
    check('locked model n_up', int((pd.to_numeric(g['sign']) > 0).sum()), m['n_up'], source='lock/diagnostic_v1/final_model_genes.tsv')
    check('locked model n_down', int((pd.to_numeric(g['sign']) < 0).sum()), m['n_down'], source='lock/diagnostic_v1/final_model_genes.tsv')
    spec = open('lock/diagnostic_v1/final_model_spec.txt').read()
    check('locked model cutoff in spec', f"cutoff = {m['cutoff']}" in spec, True, source='lock/diagnostic_v1/final_model_spec.txt')
    if a.scores:
        s = pd.read_csv(a.scores, sep='\t')
        ref = pd.read_csv(m['reference_scores'], sep='\t').set_index('sample')
        s = s.set_index('sample')
        check('scores: same 19 samples as reference', sorted(s.index) == sorted(ref.index), True, source=a.scores)
        diff = float((s.score - ref.score.reindex(s.index)).abs().max())
        check('scores: max |score - reference|', f'{diff:.2e}', f"<= {m['max_abs_score_diff']:.0e}", ok=diff <= m['max_abs_score_diff'], source=a.scores)
        same = bool((s.call_tumour.astype(int) == ref.pred_tumour.reindex(s.index).astype(int)).all())
        check('scores: calls identical to reference', same, True, source=a.scores)
        correct = int((s.call_tumour.astype(int) == ref.truth.reindex(s.index).eq('tumour').astype(int)).sum())
        check('scores: development calls correct (optimistic, not validation)', correct, m['development_correct'], source=a.scores)

    # 6. Rule 4 self-consistency (the 19-sample fit compared with itself must reproduce every candidate)
    if a.rule4:
        r4 = pd.read_csv(a.rule4, sep='\t')
        check('rule4 self-consistency: candidates evaluated', len(r4), E['final_candidates']['n'], source=a.rule4)
        check('rule4 self-consistency: reproduced', int(r4.reproduced.astype(str).eq('True').sum()), E['final_candidates']['n'], source=a.rule4)

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out.to_csv(a.out, sep='\t', index=False)
    n_fail = int(out.status.eq('FAIL').sum())
    print(out[['check', 'observed', 'expected', 'status']].to_string(index=False))
    print(f'\n{len(out)} checks: {len(out) - n_fail} PASS, {n_fail} FAIL -> {a.out}')
    sys.exit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
