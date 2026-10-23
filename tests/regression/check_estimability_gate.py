#!/usr/bin/env python3
"""Unit test of the estimability gate (upgrade H2b) in nextflow/bin/validate_analysis_input.py (check 6) and nextflow/bin/check_discover_input.py (C7).
Synthetic sample sheets (fine, nested batch, empty cell, single batch, no residual degrees of freedom) plus, read-only, the real demo and 19-sample sheets
(they must still pass). Writes only into --work (a new or empty folder). No pipeline is run.
usage: <pyqc python> -I tests/regression/check_estimability_gate.py --work <folder> [--real-sheets <00_input folder> ...]
exit 0 = all cases as expected."""
import argparse, csv, os, subprocess, sys

ap = argparse.ArgumentParser(); ap.add_argument('--work', required=True); ap.add_argument('--real-sheets', nargs='*', default=[])
a = ap.parse_args(); os.makedirs(a.work, exist_ok=True)
BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'nextflow', 'bin')
PY = sys.executable

def write_case(name, layout):
    """layout: list of (batch, group). Returns the case folder."""
    d = os.path.join(a.work, name); os.makedirs(d, exist_ok=True)
    ids = ['s%d' % i for i in range(len(layout))]
    with open(os.path.join(d, 'sheet.tsv'), 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'label', 'group', 'baseline', 'study', 'batch', 'count_source', 'patient_id'])
        for i, (b, g) in zip(ids, layout): w.writerow([i, i, g, 'TRUE' if g == 'normal' else 'FALSE', b, b, 'featurecounts', 'P' + i])
    with open(os.path.join(d, 'counts.tsv'), 'w') as f:
        f.write('gene_id\t' + '\t'.join(ids) + '\n')
        for k in range(5): f.write('G%d\t' % k + '\t'.join(str(10 + k + j) for j in range(len(ids))) + '\n')
    with open(os.path.join(d, 'symbols.tsv'), 'w') as f:
        f.write('gene_id\tgene_name\n' + ''.join('G%d\tS%d\n' % (k, k) for k in range(5)))
    return d

def run_validate(d, counts=None, sheet=None, min_samples=1, symbols=None):
    p = subprocess.run([PY, '-I', os.path.join(BIN, 'validate_analysis_input.py'), '--counts', counts or d + '/counts.tsv', '--sheet', sheet or d + '/sheet.tsv',
        '--symbols', symbols or d + '/symbols.tsv', '--route', 'featurecounts', '--min_samples', str(min_samples), '--max_unmapped_pct', '100',
        '--out_sheet', d + '/validated.tsv', '--out_json', d + '/check.json'], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr

def run_discover(d, counts=None, clinical=None, min_samples=1, symbols=None):
    cwd = os.path.join(d, 'discover_out'); os.makedirs(cwd, exist_ok=True)
    p = subprocess.run([PY, '-I', os.path.join(BIN, 'check_discover_input.py'), '--question', 'diagnostic', '--counts', counts or d + '/counts.tsv',
        '--clinical', clinical or d + '/sheet.tsv', '--symbols', symbols or d + '/symbols.tsv', '--min-samples', str(min_samples), '--max-unmapped-pct', '100'],
        capture_output=True, text=True, cwd=cwd)
    return p.returncode, p.stdout + p.stderr

N, T = 'normal', 'tumour'
cases = [  # name, layout, expect_pass, text expected in message when failing
    ('fine_two_batches', [('A', N), ('A', N), ('A', T), ('A', T), ('B', N), ('B', N), ('B', T), ('B', T)], True, ''),
    ('fine_unbalanced', [('A', N), ('A', T), ('A', T), ('A', T), ('B', N), ('B', N), ('B', N), ('B', T)], True, ''),
    ('single_batch', [('A', N), ('A', N), ('A', N), ('A', T), ('A', T), ('A', T)], True, ''),
    ('nested_batch', [('A', N), ('A', N), ('A', N), ('B', T), ('B', T), ('B', T)], False, 'batch-by-group table'),
    ('empty_cell', [('A', N), ('A', N), ('A', T), ('A', T), ('B', N), ('B', N), ('B', T), ('B', T), ('C', T), ('C', T)], False, 'batch-by-group table'),
    ('no_residual_df', [('A', N), ('A', T)], False, 'residual degrees of freedom'),
    ('minimal_ok_2x2', [('A', N), ('A', T), ('B', N), ('B', T)], True, ''),
]
bad = 0
def report(name, tool, rc, out, expect, text):
    global bad
    good = (rc == 0) == expect and (expect or (rc == 2 and text in out))
    bad += not good
    print('%s  %-18s %-9s exit %d%s' % ('PASS' if good else 'FAIL', name, tool, rc, '' if expect else '   | ' + ' '.join(out.split())[:230]))
for name, layout, expect, text in cases:
    d = write_case(name, layout)
    rc, out = run_validate(d); report(name, 'validate', rc, out, expect, text)
    rc, out = run_discover(d); report(name, 'discover', rc, out, expect, text)
for r in a.real_sheets:   # real, currently valid inputs: the 00_input folder of a DISCOVER run (read only)
    name = os.path.basename(os.path.dirname(os.path.dirname(r.rstrip('/')))) or 'real'
    d = os.path.join(a.work, 'real_' + name); os.makedirs(d, exist_ok=True)
    sym = os.path.join(d, 'symbols.tsv'); open(sym, 'w').write('gene_id\tgene_name\nx\ty\n')
    rc, out = run_validate(d, counts=r + '/counts_training.tsv', sheet=r + '/sample_sheet.discover.tsv', min_samples=3, symbols=sym); report('real:' + name, 'validate', rc, out, True, '')
    rc, out = run_discover(d, counts=r + '/counts_training.tsv', clinical=r + '/sample_sheet.discover.tsv', min_samples=3, symbols=sym); report('real:' + name, 'discover', rc, out, True, '')
print('ALL AS EXPECTED' if not bad else '%d case(s) NOT as expected' % bad); sys.exit(1 if bad else 0)
