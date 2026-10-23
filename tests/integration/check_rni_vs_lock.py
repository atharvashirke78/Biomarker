#!/usr/bin/env python3
"""Folder-entry acceptance (Addendum 6 set R19: client-style IDs RNI_001.., labels Tumor/Normal, cohorts CohortA/B) vs the LOCKED results.

usage (project root): env/prefix/pyqc/bin/python3 -I tests/integration/check_rni_vs_lock.py --discover <run> --apply <run> --baseline <pod baseline run> --out <new tsv>
Read-only. Sample IDs are mapped back with results_rc1prep_19/acceptance_format_20261016/inputs/id_map.tsv; gene IDs need no mapping.
Exit 0 only if every check is PASS."""
import argparse
import filecmp
import os
import sys

import pandas as pd

ap = argparse.ArgumentParser()
for k in ('discover', 'apply', 'baseline', 'out'):
    ap.add_argument('--' + k, required=True)
a = ap.parse_args()
if os.path.exists(a.out):
    sys.exit('ERROR: output exists (never overwritten): ' + a.out)
rows = []
add = lambda c, o, e, ok: rows.append(dict(check=c, observed=str(o)[:200], expected=str(e)[:200], status='PASS' if ok else 'FAIL'))
tsv = lambda p: pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False)
D = os.path.join(a.discover, 'pipeline', 'discover')
B = os.path.join(a.baseline, 'discover', 'out')
idm = tsv('results_rc1prep_19/acceptance_format_20261016/inputs/id_map.tsv')
new2old = dict(zip(idm.iloc[:, 1], idm.iloc[:, 0])) if idm.iloc[0, 1].startswith('RNI') else dict(zip(idm.iloc[:, 0], idm.iloc[:, 1]))

gm = tsv(os.path.join(a.discover, 'input', 'group_mapping.tsv'))
add('group mapping logged', '; '.join(f'{r.label_in_metadata}->{r.canonical_group} ({r.n_samples})' for r in gm.itertuples()), 'Tumor->tumour (9); Normal->normal (10)',
    set(zip(gm.label_in_metadata, gm.canonical_group, gm.n_samples)) == {('Tumor', 'tumour', '9'), ('Normal', 'normal', '10')})
locked = list(tsv('results_rc1prep_19/final_candidates/final_candidates_18of19.tsv').gene_id)
ct = tsv(os.path.join(D, '02_candidates', 'candidate_table.tsv'))
new = list(ct.gene_id[ct.all_rules_pass == 'TRUE'])
add('candidates: genes and order = final_candidates_18of19.tsv', len(new), f'{len(locked)} same order', new == locked)
for f in ('01_de/deseq2_results_all.tsv', '01_de/deseq2_salmon_results_all.tsv', '01_de/deseq2_summary.json',
          '02_candidates/funnel_counts.tsv', '03_model/model_genes.tsv', '03_model/model_spec.txt', '03_model/inner_fold_cutoffs.tsv'):
    same = filecmp.cmp(os.path.join(D, f), os.path.join(B, f), shallow=False)
    add(f'byte-identical to pod baseline (SRR IDs): {f}', same, True, same)
# candidate_table.tsv names the left-out samples a gene failed on (sample IDs): compare after mapping RNI -> SRR, nothing else changed
txt = open(os.path.join(D, '02_candidates', 'candidate_table.tsv')).read()
for n, o in sorted(new2old.items(), key=lambda kv: -len(kv[0])):
    txt = txt.replace(n, o)
same = txt == open(os.path.join(B, '02_candidates', 'candidate_table.tsv')).read()
add('identical to pod baseline after mapping sample IDs back: 02_candidates/candidate_table.tsv', same, True, same)
lg = tsv('lock/diagnostic_v1/final_model_genes.tsv').set_index('gene_id')
mg = tsv(os.path.join(D, '03_model', 'model_genes.tsv')).set_index('gene_id')
add('model genes, signs, coefficients = lock', sorted(mg.index), sorted(lg.index),
    set(mg.index) == set(lg.index) and all(mg.loc[g, 'sign'] == lg.loc[g, 'sign'] and mg.loc[g, 'lasso_coef_lambda1se'] == lg.loc[g, 'lasso_coef_lambda1se'] for g in lg.index))
cut = lambda p: [l.split()[2] for l in open(p) if l.startswith('cutoff = ')]
add('cutoff = lock', cut(os.path.join(D, '03_model', 'model_spec.txt')), cut('lock/diagnostic_v1/final_model_spec.txt'),
    cut(os.path.join(D, '03_model', 'model_spec.txt')) == cut('lock/diagnostic_v1/final_model_spec.txt'))
ref = pd.read_csv('results_rc1prep_19/lock_proposal_v2/development_scores.tsv', sep='\t').set_index('sample')
for name, p, col, call in (('DISCOVER development scores', os.path.join(D, '03_model', 'development_scores.tsv'), 'score', 'pred_positive'),
                           ('APPLY locked-model scores', os.path.join(a.apply, 'pipeline', 'apply', 'scores', 'scores.tsv'), 'score', 'call_tumour')):
    s = pd.read_csv(p, sep='\t')
    s['old'] = s['sample'].map(new2old)
    s = s.set_index('old')
    d = (s[col] - ref.score.reindex(s.index)).abs().max()
    add(f'{name} (RNI -> SRR) = locked development scores, max abs diff', d, 0, set(s.index) == set(ref.index) and d == 0)
    add(f'{name}: calls = locked calls', bool((s[call] == ref.pred_tumour.reindex(s.index)).all()), True, bool((s[call] == ref.pred_tumour.reindex(s.index)).all()))
out = pd.DataFrame(rows)
out.to_csv(a.out, sep='\t', index=False)
print(out[['status', 'check', 'observed']].to_string(index=False, max_colwidth=80))
print(f"\n{(out.status == 'PASS').sum()} of {len(out)} PASS")
sys.exit(0 if (out.status == 'PASS').all() else 1)
