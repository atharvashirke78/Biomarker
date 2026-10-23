#!/usr/bin/env python3
"""Two-level regression, level 1 (the "second check script"): decision-level comparison with the LOCKED results (read-only).

Independent of the pipeline code: it reads the lock files and the recorded locked-era results directly and recomputes the
locked-model score itself (Python, from the counts matrix and lock/diagnostic_v1), instead of trusting the pipeline's numbers.

usage (project root): env/prefix/pyqc/bin/python3 -I tests/regression/check_decisions.py --run <dir from run_regression.sh> --out <new tsv>
References (never written):
  lock/diagnostic_v1/{sha256sums_final.txt, final_model_genes.tsv, final_model_spec.txt}
  results_rc1prep_19/final_candidates/final_candidates_18of19.tsv          the 1,390 locked candidates, in order
  results_rc1prep_19/discover_regression_20261015/out/                     locked-era DISCOVER outputs (rule flags, DE significant sets)
  results_rc1prep_19/lock_proposal_v2/development_scores.tsv               development scores of the locked model
  results/04_counts/counts_featurecounts.tsv                               counts for the independent score recomputation
Checks: lock hashes; candidate list (genes, order); every rule flag; DE significant sets (both routes); model genes, signs,
coefficients (as printed in the lock, 6 significant digits), cutoff; DISCOVER development scores and calls; locked-model scores
and calls from run_diagnostic.sh; independent recomputation of the locked scores; rule 4 self-consistency 1,390 of 1,390.
Exit 0 only if every check is PASS. Nothing is rounded or changed to make a check pass beyond what is stated per check."""
import argparse
import hashlib
import json
import math
import os
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--run', required=True)
ap.add_argument('--out', required=True)
a = ap.parse_args()
if os.path.exists(a.out):
    sys.exit('ERROR: output exists (never overwritten): ' + a.out)
R = 'results_rc1prep_19'
LOCK = 'lock/diagnostic_v1'
D = os.path.join(a.run, 'discover', 'out')
G = os.path.join(a.run, 'diagnostic')
OLD = f'{R}/discover_regression_20261015/out'
rows = []


def add(check, observed, expected, ok):
    rows.append(dict(check=check, observed=str(observed)[:300], expected=str(expected)[:300], status='PASS' if ok else 'FAIL'))


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def tsv(p):
    return pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False)


# 1. lock integrity
bad = []
for line in open(f'{LOCK}/sha256sums_final.txt'):
    h, p = line.split()
    if sha(p) != h:
        bad.append(p)
add('lock hashes (sha256sums_final.txt)', f'{7 - len(bad)} of 7 OK', '7 of 7 OK', not bad)

# 2. candidates
fc = tsv(f'{R}/final_candidates/final_candidates_18of19.tsv')
locked_ids = list(fc.gene_id)
ct = tsv(f'{D}/02_candidates/candidate_table.tsv')
new_ids = list(ct.gene_id[ct.all_rules_pass == 'TRUE'])
add('candidates passing all rules: genes and order = final_candidates_18of19.tsv', f'{len(new_ids)}', f'{len(locked_ids)} (same order)', new_ids == locked_ids)
ot = tsv(f'{OLD}/02_candidates/candidate_table.tsv').set_index('gene_id')
nt = ct.set_index('gene_id')
add('candidate table: same gene universe as locked-era run', len(nt), len(ot), set(nt.index) == set(ot.index))
for c in ('rule1_fc', 'rule1_salmon', 'direction_fc', 'direction_salmon', 'rule1_pass', 'rule2_pass', 'baseMean_ge_50_both', 'rule3a_pass', 'loo_n_pass', 'loo_n_fail', 'rule3b_pass', 'all_rules_pass'):
    if c in ot.columns:
        k = int((ot[c] != nt.loc[ot.index, c]).sum())
        add(f'rule flag {c}: genes differing from locked-era run', k, 0, k == 0)

# 3. DE significant sets
for f in ('deseq2_significant_up.tsv', 'deseq2_significant_down.tsv', 'deseq2_salmon_significant_up.tsv', 'deseq2_salmon_significant_down.tsv'):
    so, sn = set(tsv(f'{OLD}/01_de/{f}').gene_id), set(tsv(f'{D}/01_de/{f}').gene_id)
    add(f'DE significant set {f}', f'{len(sn)} (symmetric diff {len(so ^ sn)})', len(so), so == sn)

# 4. model vs lock
lg = tsv(f'{LOCK}/final_model_genes.tsv').set_index('gene_id')
mg = tsv(f'{D}/03_model/model_genes.tsv').set_index('gene_id')
add('model genes = lock', sorted(mg.index), sorted(lg.index), set(mg.index) == set(lg.index))
if set(mg.index) == set(lg.index):
    s_ok = all(int(mg.loc[g, 'sign']) == int(lg.loc[g, 'sign']) for g in lg.index)
    add('model signs = lock', 'all equal' if s_ok else 'differ', 'all equal', s_ok)
    c_ok = all(mg.loc[g, 'lasso_coef_lambda1se'] == lg.loc[g, 'lasso_coef_lambda1se'] for g in lg.index)
    add('LASSO coefficients = lock (as printed, 6 significant digits)', 'all equal' if c_ok else dict(mg.lasso_coef_lambda1se), 'all equal', c_ok)
cut_lock = [l.split()[2] for l in open(f'{LOCK}/final_model_spec.txt') if l.startswith('cutoff = ')]
cut_new = [l.split()[2] for l in open(f'{D}/03_model/model_spec.txt') if l.startswith('cutoff = ')]
add('cutoff = lock', cut_new, cut_lock, len(cut_lock) == 1 and cut_new == cut_lock)

# 5. scores
ref = pd.read_csv(f'{R}/lock_proposal_v2/development_scores.tsv', sep='\t').set_index('sample')
ds = pd.read_csv(f'{D}/03_model/development_scores.tsv', sep='\t').set_index('sample')
d1 = (ds.score - ref.score.reindex(ds.index)).abs().max()
add('DISCOVER development scores = locked development scores (max abs diff)', d1, 0, set(ds.index) == set(ref.index) and d1 == 0)
add('DISCOVER development calls = locked calls', bool((ds.pred_positive == ref.pred_tumour.reindex(ds.index)).all()), True, bool((ds.pred_positive == ref.pred_tumour.reindex(ds.index)).all()))
ls = pd.read_csv(f'{G}/locked_model_development/scores.tsv', sep='\t').set_index('sample')
d2 = (ls.score - ref.score.reindex(ls.index)).abs().max()
add('apply_locked_model.R scores = locked development scores (max abs diff)', d2, 0, set(ls.index) == set(ref.index) and d2 == 0)
add('apply_locked_model.R calls = locked calls', bool((ls.call_tumour == ref.pred_tumour.reindex(ls.index)).all()), True, bool((ls.call_tumour == ref.pred_tumour.reindex(ls.index)).all()))

# 6. independent recomputation of the locked score (does not use any pipeline script)
cut = float(cut_lock[0])
cn = pd.read_csv('results/04_counts/counts_featurecounts.tsv', sep='\t', index_col=0)
tot = cn.sum(axis=0)
up = [g for g in lg.index if int(lg.loc[g, 'sign']) > 0]
dn = [g for g in lg.index if int(lg.loc[g, 'sign']) < 0]
lc = (cn.loc[up + dn] / tot * 1e6 + 1).apply(lambda col: col.map(math.log2))
mine = lc.loc[up].mean() - lc.loc[dn].mean()
d3 = (mine - ls.score.reindex(mine.index)).abs().max()
add('independent Python recomputation vs apply_locked_model.R (max abs diff; floating point only)', f'{d3:.3g}', '< 1e-9', d3 < 1e-9)
calls = (mine > cut).astype(int)
add('independent recomputation: calls = locked calls', bool((calls == ref.pred_tumour.reindex(calls.index)).all()), True, bool((calls == ref.pred_tumour.reindex(calls.index)).all()))

# 7. rule 4 self-consistency
r4 = json.load(open(f'{G}/rule4_self_consistency/rule4_summary.json'))
add('rule 4 self-consistency (development fit)', f"{r4['n_reproduced']} of {r4['n_candidates']}", '1390 of 1390', r4['n_reproduced'] == 1390 and r4['n_candidates'] == 1390)

out = pd.DataFrame(rows)
out.to_csv(a.out, sep='\t', index=False)
print(out[['status', 'check', 'observed']].to_string(index=False, max_colwidth=90))
print(f"\n{(out.status == 'PASS').sum()} of {len(out)} PASS")
sys.exit(0 if (out.status == 'PASS').all() else 1)
