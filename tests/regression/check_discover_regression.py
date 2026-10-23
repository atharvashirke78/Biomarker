#!/usr/bin/env python3
"""Regression: -entry DISCOVER on the 19 development samples vs the existing (locked) results. Read-only on everything it compares.

usage: env/prefix/pyqc/bin/python3 -I tests/regression/check_discover_regression.py --out-dir <DISCOVER outdir> --summary <new summary.tsv>
Each row: item, comparison, observed, expected, status (IDENTICAL / DIFFERENT / PASS / FAIL). Exit 0 only if nothing is DIFFERENT or FAIL.
Nothing is changed to make a comparison match: a difference is reported as it is.
"""
import argparse
import filecmp
import os
import re
import sys

import pandas as pd

rows = []


def add(item, comparison, observed, expected, ok, kind='byte'):
    status = ('IDENTICAL' if ok else 'DIFFERENT') if kind == 'byte' else ('PASS' if ok else 'FAIL')
    rows.append(dict(item=item, comparison=comparison, observed=observed, expected=expected, status=status))


def same(a, b, item):
    if not (os.path.exists(a) and os.path.exists(b)):
        add(item, f'{a} vs {b}', 'missing' if not os.path.exists(a) else 'present', 'present', False)
        return False
    ok = filecmp.cmp(a, b, shallow=False)
    add(item, f'byte compare with {b}', 'identical' if ok else 'differs', 'identical', ok)
    return ok


ap = argparse.ArgumentParser()
ap.add_argument('--out-dir', required=True)
ap.add_argument('--summary', required=True)
a = ap.parse_args()
O = a.out_dir
R = 'results_rc1prep_19'
LOCK = 'lock/diagnostic_v1'

# ---- DE (both routes)
for f in ('deseq2_results_all.tsv', 'deseq2_summary.json', 'deseq2_salmon_results_all.tsv', 'deseq2_salmon_summary.json'):
    same(f'{O}/01_de/{f}', f'{R}/{f}', f'DE {f}')

# ---- candidate rules
C = f'{O}/02_candidates'
for f in ('candidate_table.tsv', 'loo_pass_matrix.tsv', 'loo_failures_by_dropped_sample.tsv', 'reproduction_check_19_sample_fit.tsv'):
    same(f'{C}/{f}', f'{R}/rule3b_rerun18_cfg/{f}', f'candidates {f}')
same(f'{C}/funnel_counts.tsv', f'{R}/diagnostic_final_20261015/regression/candidate_selection_rerun/funnel_counts.tsv', 'candidates funnel_counts.tsv (corrected label "18 of 19")')
ct = pd.read_csv(f'{C}/candidate_table.tsv', sep='\t')
fin = ct[ct.all_rules_pass.astype(str) == 'True']
locked_final = pd.read_csv(f'{R}/final_candidates/final_candidates_18of19.tsv', sep='\t')
add('candidates n', 'all_rules_pass rows', len(fin), 1390, len(fin) == 1390, 'num')
add('candidates up / down', 'direction_fc', f"{int((fin.direction_fc == 'up').sum())} / {int((fin.direction_fc == 'down').sum())}", '510 / 880',
    (fin.direction_fc == 'up').sum() == 510 and (fin.direction_fc == 'down').sum() == 880, 'num')
add('candidates protein-coding', 'gene_biotype', int((fin.gene_biotype == 'protein_coding').sum()), 1353, (fin.gene_biotype == 'protein_coding').sum() == 1353, 'num')
add('candidates = final_candidates_18of19.tsv', 'same gene IDs in the same order', list(fin.gene_id) == list(locked_final.gene_id), True, list(fin.gene_id) == list(locked_final.gene_id), 'num')

# ---- model vs lock
M = f'{O}/03_model'
mg = pd.read_csv(f'{M}/model_genes.tsv', sep='\t')
lg = pd.read_csv(f'{LOCK}/final_model_genes.tsv', sep='\t')
add('model genes', 'gene_id set vs lock', len(mg), 9, set(mg.gene_id) == set(lg.gene_id), 'num')
sg = dict(zip(mg.gene_id, mg.set)) == dict(zip(lg.gene_id, lg.set))
add('model signs', 'up/down per gene vs lock', 'same' if sg else 'differs', 'same', sg, 'num')
coef = pd.Series(dict(zip(mg.gene_id, mg.lasso_coef_lambda1se))); lcoef = pd.Series(dict(zip(lg.gene_id, lg.lasso_coef_lambda1se)))
d = float((coef - lcoef.reindex(coef.index)).abs().max()) if set(coef.index) == set(lcoef.index) else float('nan')
add('model coefficients', 'max |coef - lock coef| (both signif 6)', f'{d:.3g}', '0', d == 0, 'num')
spec = open(f'{M}/model_spec.txt').read(); lspec = open(f'{LOCK}/final_model_spec.txt').read()
cut = re.search(r'^cutoff = (\S+)', spec, re.M).group(1); lcut = re.search(r'^cutoff = (\S+)', lspec, re.M).group(1)
add('cutoff', 'model_spec "cutoff =" vs lock (%.10g text)', cut, lcut, cut == lcut, 'num')
add('cutoff numeric', '|cutoff - 1.352325902|', f'{abs(float(cut) - 1.352325902):.2e}', '<= 1e-9', abs(float(cut) - 1.352325902) <= 1e-9, 'num')
fo = re.search(r'Youden cutoffs: ([^)]+)\)', spec).group(1); lfo = re.search(r'Youden cutoffs: ([^)]+)\)', lspec).group(1)
add('inner-fold cutoffs', 'as printed in the spec', fo, lfo, fo == lfo, 'num')
lam = re.search(r'lambda\.1se = (\S+)', spec).group(1); llam = re.search(r'lambda\.1se = (\S+)', lspec).group(1)
add('lambda.1se', 'as printed in the spec (%.8g)', lam, llam, lam == llam, 'num')
sc = pd.read_csv(f'{M}/development_scores.tsv', sep='\t').set_index('sample')
ls = pd.read_csv(f'{R}/lock_proposal_v2/development_scores.tsv', sep='\t').set_index('sample')
ds = float((sc.score - ls.score.reindex(sc.index)).abs().max())
add('development scores', 'max |score - lock_proposal_v2 score|', f'{ds:.2e}', '<= 1e-10', ds <= 1e-10 and set(sc.index) == set(ls.index), 'num')
calls = bool((sc.pred_positive == ls.pred_tumour.reindex(sc.index)).all())
add('development calls', 'identical to lock_proposal_v2', calls, True, calls, 'num')
pf = pd.read_csv(f'{M}/development_performance.tsv', sep='\t').iloc[0]
add('development performance (OPTIMISTIC)', 'AUC, TP/FN/TN/FP', f"{pf.AUC:g}, {pf.TP}/{pf.FN}/{pf.TN}/{pf.FP}", '1, 9/0/10/0', (pf.AUC, pf.TP, pf.FN, pf.TN, pf.FP) == (1, 9, 0, 10, 0), 'num')

out = pd.DataFrame(rows)
out.to_csv(a.summary, sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
print(out[['item', 'observed', 'expected', 'status']].to_string(index=False))
bad = out.status.isin(['DIFFERENT', 'FAIL']).sum()
print(f'\n{len(out)} comparisons: {len(out) - bad} identical/pass, {bad} different/fail -> {a.summary}')
sys.exit(1 if bad else 0)
