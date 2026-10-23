#!/usr/bin/env python3
"""Rule 4 (LOCKED, config/candidate_rules.yaml rule4_independent_reproduction): independent reproduction of the 1,390 candidates.

A candidate gene is "reproduced" if, in the independent samples, it has the same direction as in the 19 development samples
(sign of the shrunk log2FC, featureCounts route, in both) AND padj < alpha (de_padj, 0.05) in a DESeq2 run that uses the same
model as the 19-sample run: ~ batch + group if more than one study is present, ~ group otherwise, same pre-filter rule
(>= 10 counts in >= n samples, n = smallest group). That DESeq2 run is produced by the existing nextflow/bin/deseq2.R
(e.g. nextflow -entry ANALYSIS on the QUANT contract files); this script only compares. Reported: number and percentage of
the candidates reproduced. No pass or fail threshold (pass_threshold: none).

usage (project root, pyqc env):
  env/prefix/pyqc/bin/python3 -I scripts/rule4_reproduction.py --results <deseq2_results_all.tsv> --summary <deseq2_summary.json>
      --candidates results_rc1prep_19/final_candidates/final_candidates_18of19.tsv --label <independent|self-consistency> --out-dir <new dir>
--label self-consistency is for tests with the development fit itself (every candidate must be reproduced); it is NOT validation.
"""
import argparse
import hashlib
import json
import os
import re
import sys

import pandas as pd


def die(msg):
    print('ERROR: ' + msg, file=sys.stderr)
    sys.exit(2)


def read_key(path, key):
    for line in open(path):
        m = re.match(r'\s*' + re.escape(key) + r'\s*[:=]\s*([^\s#/]+)', line)
        if m:
            return m.group(1).strip('"\'')
    die(f'{key} not found in {path}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', required=True)
    ap.add_argument('--summary', required=True)
    ap.add_argument('--candidates', required=True)
    ap.add_argument('--label', required=True, choices=['independent', 'self-consistency'])
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    if os.path.exists(a.out_dir):
        die(f'output folder exists (never overwritten): {a.out_dir}')

    rules = open('config/candidate_rules.yaml').read()
    m = re.search(r'^rule4_independent_reproduction:.*?\n\s+status:\s*(\S+)', rules, re.M | re.S)
    if not m or m.group(1) != 'LOCKED':
        die('rule4_independent_reproduction with status LOCKED not found in config/candidate_rules.yaml')
    alpha = float(read_key('conf/thresholds.config', 'de_padj'))

    summ = json.load(open(a.summary))
    design, nb = summ.get('design', ''), None
    res = pd.read_csv(a.results, sep='\t')
    for col in ('gene_id', 'padj', 'log2FoldChange_shrunk'):
        if col not in res.columns:
            die(f'results table lacks column {col}')
    if summ.get('padj_cut') is not None and abs(float(summ['padj_cut']) - alpha) > 1e-12:
        die(f"DESeq2 run used padj_cut {summ['padj_cut']}, config de_padj is {alpha}")
    if design not in ('~ batch + group', '~ group'):
        die(f'unexpected DESeq2 design in summary: {design!r}')

    cand = pd.read_csv(a.candidates, sep='\t')
    if 'log2FC_shrunk_fc' not in cand.columns:
        die('candidates table lacks log2FC_shrunk_fc')
    r = res.set_index('gene_id')
    t = cand[['gene_id', 'symbol', 'direction_fc', 'log2FC_shrunk_fc']].copy()
    t['in_independent_table'] = t.gene_id.isin(r.index)
    t['indep_log2FC_shrunk'] = t.gene_id.map(r.log2FoldChange_shrunk)
    t['indep_padj'] = t.gene_id.map(r.padj)
    same = (t.indep_log2FC_shrunk * t.log2FC_shrunk_fc) > 0
    t['same_direction'] = same
    t['reproduced'] = same & t.indep_padj.lt(alpha)
    t['status'] = 'reproduced'
    t.loc[~t.in_independent_table, 'status'] = 'not tested (removed by pre-filter or absent)'
    t.loc[t.in_independent_table & t.indep_padj.isna(), 'status'] = 'not tested (padj NA)'
    t.loc[t.in_independent_table & t.indep_padj.notna() & ~t.reproduced & ~same, 'status'] = 'opposite direction'
    t.loc[t.in_independent_table & t.indep_padj.notna() & ~t.reproduced & same, 'status'] = f'same direction, padj >= {alpha}'

    n, k = len(t), int(t.reproduced.sum())
    note = ('SELF-CONSISTENCY TEST on the development fit: not validation' if a.label == 'self-consistency'
            else 'independent samples, locked rule 4, reported only (no pass threshold)')
    os.makedirs(a.out_dir)
    t.to_csv(os.path.join(a.out_dir, 'rule4_per_gene.tsv'), sep='\t', index=False)
    s = t.status.value_counts().rename_axis('status').reset_index(name='n')
    s.to_csv(os.path.join(a.out_dir, 'rule4_status_counts.tsv'), sep='\t', index=False)
    sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
    info = {'label': a.label, 'note': note, 'candidates': a.candidates, 'candidates_sha256': sha(a.candidates),
            'results': a.results, 'results_sha256': sha(a.results), 'design': design, 'n_samples': summ.get('n_samples'),
            'alpha': alpha, 'n_candidates': n, 'n_reproduced': k, 'pct_reproduced': round(100 * k / n, 2),
            'by_direction': {d: {'n': int((t.direction_fc == d).sum()), 'reproduced': int((t.reproduced & t.direction_fc.eq(d)).sum())} for d in ('up', 'down')},
            'pass_threshold': 'none (reported only)', 'command': ' '.join(sys.argv)}
    json.dump(info, open(os.path.join(a.out_dir, 'rule4_summary.json'), 'w'), indent=2)
    print(f'rule 4 [{note}]: {k} of {n} candidates reproduced ({100 * k / n:.1f} %); design {design}; n_samples {summ.get("n_samples")}')
    print(s.to_string(index=False))


if __name__ == '__main__':
    main()
