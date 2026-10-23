#!/usr/bin/env python3
"""Compare two pipeline output folders (new vs baseline) stage by stage, read-only. Prints one line per check and exits 1 if a check that must be exact differs.

usage: <pyqc python> tests/regression/compare_runs.py --new <output folder> --base <output folder> [--tsv result.tsv]
Checks: featureCounts matrix (sha256), sample QC, DESeq2 featureCounts route, leave-one-out matrix, candidate columns of the featureCounts route, funnel, final candidates
(must be EXACT), and the Salmon route (expected to differ a little: Salmon is not exactly repeatable; reported with sizes, never a failure).
Tolerance for the exact checks is 0 (cells are read with round-trip float parsing); the Salmon checks only report."""
import argparse
import hashlib
import os
import sys

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--new', required=True)
ap.add_argument('--base', required=True)
ap.add_argument('--tsv')
a = ap.parse_args()
N, B = a.new, a.base
rows = []


def add(check, status, detail):
    rows.append((check, status, detail))
    print(f'{status:10s} {check:52s} {detail}')


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()[:16] if os.path.exists(p) else None


def table(base, rel, key=None):
    p = os.path.join(base, rel)
    if not os.path.exists(p):
        return None
    t = pd.read_csv(p, sep='\t', float_precision='round_trip')
    return t.set_index(key) if key else t


def same_table(check, rel, key, cols=None, must=True):
    n, b = table(N, rel, key), table(B, rel, key)
    if n is None or b is None:
        return add(check, 'MISSING', f'{rel}: ' + ('new' if n is None else 'baseline') + ' file not found')
    if cols:
        n, b = n[cols], b[cols]
    if n.shape != b.shape or list(n.columns) != list(b.columns) or not n.index.equals(b.index):
        only = sorted(set(n.index) ^ set(b.index))[:5]
        return add(check, 'DIFFERENT' if must else 'differs', f'shape new {n.shape} base {b.shape}; rows only in one: {only}')
    bad = []
    for c in n.columns:
        x, y = n[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y) and not pd.api.types.is_bool_dtype(x):
            d = np.abs(x.values.astype(float) - y.values.astype(float))
            nan_mismatch = (np.isnan(x.values.astype(float)) != np.isnan(y.values.astype(float))).sum()
            if nan_mismatch or (np.nanmax(d) > 0 if np.isfinite(d).any() else False):
                bad.append(f'{c} (max |diff| {np.nanmax(d) if np.isfinite(d).any() else 0:.3g})')
        elif not (x.fillna('').astype(str) == y.fillna('').astype(str)).all():
            bad.append(c + ' (text)')
    if bad:
        return add(check, 'DIFFERENT' if must else 'differs', 'columns differ: ' + ', '.join(bad[:6]))
    add(check, 'EXACT', f'{n.shape[0]} rows x {n.shape[1]} columns, every cell equal')


a_, b_ = sha(f'{N}/04_counts/counts_featurecounts.tsv'), sha(f'{B}/04_counts/counts_featurecounts.tsv')
add('featureCounts matrix (sha256)', 'EXACT' if a_ and a_ == b_ else ('MISSING' if not (a_ and b_) else 'DIFFERENT'), f'new {a_} base {b_}')
same_table('sample QC table', '04_counts/sample_qc.tsv', 'srr')
same_table('DESeq2 featureCounts (all QC-pass samples)', '06_dge/deseq2_results_all.tsv', 'gene_id')
same_table('DESeq2 featureCounts (discovery)', '09_discover/01_de/deseq2_results_all.tsv', 'gene_id')
same_table('leave-one-out pass matrix', '09_discover/02_candidates/loo_pass_matrix.tsv', 'gene_id')
fc_cols = ['baseMean_fc', 'log2FC_shrunk_fc', 'padj_fc', 'rule1_fc', 'loo_n_pass']
nct, bct = table(N, '09_discover/02_candidates/candidate_table.tsv', 'gene_id'), table(B, '09_discover/02_candidates/candidate_table.tsv', 'gene_id')
if nct is not None and bct is not None:
    common = nct.index.intersection(bct.index)
    d = [c for c in fc_cols if c in nct and c in bct and not (nct.loc[common, c].astype(str) == bct.loc[common, c].astype(str)).all()]
    add('candidate table: featureCounts-route columns', 'EXACT' if not d else 'DIFFERENT', f'{len(common)} shared genes' + (f'; differ: {d}' if d else ''))
    fa, fb = set(nct.index[nct.all_rules_pass.astype(str) == 'TRUE']), set(bct.index[bct.all_rules_pass.astype(str) == 'TRUE'])
    add('final candidate genes (all rules passed)', 'EXACT' if fa == fb else 'DIFFERENT', f'new {len(fa)} base {len(fb)}; only new {sorted(fa - fb)[:5]}; only base {sorted(fb - fa)[:5]}')
    ra, rb = set(nct.index), set(bct.index)
    add('candidate-table rows (rule 1 on either route)', 'EXACT' if ra == rb else 'differs', f'new {len(ra)} base {len(rb)}; only new {sorted(ra - rb)[:3]}; only base {sorted(rb - ra)[:3]} (depends on the Salmon route)')
fn, fb_ = table(N, '09_discover/02_candidates/funnel_counts.tsv'), table(B, '09_discover/02_candidates/funnel_counts.tsv')
if fn is not None and fb_ is not None:
    m = pd.concat([fb_.iloc[:, :2].add_prefix('base_').reset_index(drop=True), fn.iloc[:, :2].add_prefix('new_').reset_index(drop=True)], axis=1)
    diff_rows = m[(m.iloc[:, 0] != m.iloc[:, 2])]
    add('candidate funnel counts', 'EXACT' if diff_rows.empty else 'differs', 'all stages equal' if diff_rows.empty else 'stages that differ: ' + '; '.join(f"{fn.iloc[i, 0] if fn.shape[1] > 2 else i}: base {m.iloc[i, 0]} new {m.iloc[i, 2]}" for i in diff_rows.index))
# Salmon route: report only
sn, sb = table(N, '04_counts/counts_salmon_tximport.tsv'), table(B, '04_counts/counts_salmon_tximport.tsv')
if sn is not None and sb is not None:
    sn, sb = sn.set_index(sn.columns[0]), sb.set_index(sb.columns[0])
    if sn.shape == sb.shape:
        l1 = ((sn - sb).abs().sum() / sb.sum() * 100)
        many = ((sn - sb).abs() > 1).sum()
        add('Salmon gene counts (expected to differ a little)', 'REPORT', f'L1 difference per sample {l1.min():.3f}-{l1.max():.3f}% of counts; genes differing by >1 read per sample {int(many.min())}-{int(many.max())}')
    else:
        add('Salmon gene counts', 'REPORT', f'shape new {sn.shape} base {sb.shape}')
dn, db_ = table(N, '09_discover/01_de/deseq2_salmon_results_all.tsv', 'gene_id'), table(B, '09_discover/01_de/deseq2_salmon_results_all.tsv', 'gene_id')
if dn is not None and db_ is not None:
    c = dn.index.intersection(db_.index)
    d = (dn.loc[c, 'log2FoldChange_shrunk'] - db_.loc[c, 'log2FoldChange_shrunk']).abs()
    add('DESeq2 Salmon route (expected to differ a little)', 'REPORT', f'{len(c)} shared genes; max |shrunk log2FC diff| {d.max():.3g}; genes only in one table: {len(set(dn.index) ^ set(db_.index))}')
must_fail = [r for r in rows if r[1] in ('DIFFERENT', 'MISSING')]
print(f'\n{len([r for r in rows if r[1] == "EXACT"])} exact, {len(must_fail)} must-be-exact checks failed or missing, {len([r for r in rows if r[1] in ("differs", "REPORT")])} informational')
if a.tsv:
    pd.DataFrame(rows, columns=['check', 'status', 'detail']).to_csv(a.tsv, sep='\t', index=False)
sys.exit(1 if must_fail else 0)
