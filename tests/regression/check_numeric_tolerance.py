#!/usr/bin/env python3
"""Acceptance check of a regression run against a CONTROL run of the frozen code on the same pod (read-only).

ACCEPTANCE RULE (owner decision 2026-10-21, goes into the freeze record): (1) every decision is identical (tests/regression/check_decisions.py:
all checks PASS on the run), and (2) every numeric value of every compared table is within a relative difference of 1e-5 of the control
(|a - b| / max(|a|, |b|) <= 1e-5; both 0 or both NaN count as equal), with identical column names, row counts, gene/sample order and identical text cells.

usage: check_numeric_tolerance.py --control <run dir> --run <run dir> --out <new tsv> [--rtol 1e-5]
Both dirs come from tests/regression/run_regression*.sh. Compared: .tsv/.csv (cell by cell), .json (numbers within the tolerance, everything else equal),
.txt/.md of results (text equal). Not compared and listed as SKIP: Nextflow work/ and .nextflow, logs, traces, reports, provenance and manifest files
(they contain run names, paths and sha256 of the files compared here), PNG/binary files, stub_default/ (a stub run has no numbers).
Also not compared: files written by the check scripts themselves (check_*.tsv / check_*.out). Text cells are compared after removing the absolute install
prefix (/workspace/Biomarker_pipeline/) and the two run-folder prefixes, because the branch copy of the scripts passes absolute paths where the original passes
relative ones; only that prefix is ignored, never a number or a name.
Exit 0 only if no file is OUT_OF_TOL / STRUCT_DIFF / MISSING / EXTRA."""
import argparse
import csv
import json
import math
import os
import re
import sys

ap = argparse.ArgumentParser()
ap.add_argument('--control', required=True)
ap.add_argument('--run', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--rtol', type=float, default=1e-5)
ap.add_argument('--allow-extra', nargs='*', default=['robustness_refits.tsv', 'development_sens_spec.tsv', 'evaluation_sens_spec.tsv'], help='file names that this branch adds and the control does not have (new outputs, not changes)')
a = ap.parse_args()
if os.path.exists(a.out):
    sys.exit('ERROR: output exists (never overwritten): ' + a.out)
C, R = os.path.realpath(a.control), os.path.realpath(a.run)
SKIP_DIR = re.compile(r'(^|/)(work|\.nextflow|stub_default)(/|$)')
SKIP_FILE = re.compile(r'(\.log(\.\d+)?$|^\.nextflow|trace.*\.tsv$|report.*\.html$|timeline.*\.html$|^exit_code\.txt$|\.console$|run_info.*|\.provenance\.txt$|^manifest.*\.json$|'
                       r'^check_.*\.(tsv|out)$|^MODEL_MANIFEST\.tsv$|^RESULTS_INDEX.*\.md$|^DISCOVER_SUMMARY\.md$|^sessionInfo\.txt$|^tool_versions\.txt$|\.png$|\.rds$|\.bam$|\.bai$|exit_code\.txt$)')


def files(root):
    out = {}
    for dp, dn, fs in os.walk(root):
        for f in fs:
            rel = os.path.relpath(os.path.join(dp, f), root)
            out[rel] = os.path.join(dp, f)
    return out


def rel(x, y):
    if x == y or (math.isnan(x) and math.isnan(y)):
        return 0.0
    if math.isnan(x) or math.isnan(y):
        return math.inf
    d = max(abs(x), abs(y))
    return abs(x - y) / d if d else 0.0


PREFIXES = [C + '/', R + '/', '/workspace/Biomarker_pipeline/']      # longest first


def norm(t):
    for pf in PREFIXES:
        t = t.replace(pf, '')
    return t


def num(s):
    try:
        return float(s)
    except ValueError:
        return None


def cmp_table(p, q, delim):
    x = list(csv.reader(open(p, newline='', errors='replace'), delimiter=delim))
    y = list(csv.reader(open(q, newline='', errors='replace'), delimiter=delim))
    if len(x) != len(y):
        return 'STRUCT_DIFF', f'rows {len(x)} vs {len(y)}', 0.0
    mx = 0.0
    for i, (r1, r2) in enumerate(zip(x, y)):
        if len(r1) != len(r2):
            return 'STRUCT_DIFF', f'line {i + 1}: {len(r1)} vs {len(r2)} fields', 0.0
        for j, (u, v) in enumerate(zip(r1, r2)):
            if u == v or norm(u) == norm(v):
                continue
            a_, b_ = num(u), num(v)
            if a_ is None or b_ is None:
                return 'STRUCT_DIFF', f'line {i + 1} field {j + 1}: text differs ({u[:40]!r} vs {v[:40]!r})', mx
            d = rel(a_, b_)
            mx = max(mx, d)
            if d > a.rtol:
                return 'OUT_OF_TOL', f'line {i + 1} field {j + 1}: {u} vs {v} (rel {d:.3g})', mx
    return ('IDENTICAL' if mx == 0 else 'WITHIN_TOL'), '', mx


def cmp_json(x, y, path=''):
    mx = 0.0
    if isinstance(x, dict) and isinstance(y, dict):
        if set(x) != set(y):
            return None, f'{path}: keys differ'
        for k in x:
            m, why = cmp_json(x[k], y[k], f'{path}/{k}')
            if m is None:
                return None, why
            mx = max(mx, m)
        return mx, ''
    if isinstance(x, list) and isinstance(y, list):
        if len(x) != len(y):
            return None, f'{path}: list length differs'
        for i, (u, v) in enumerate(zip(x, y)):
            m, why = cmp_json(u, v, f'{path}[{i}]')
            if m is None:
                return None, why
            mx = max(mx, m)
        return mx, ''
    if isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
        d = rel(float(x), float(y))
        return (d, '') if d <= a.rtol else (None, f'{path}: {x} vs {y} (rel {d:.3g})')
    return (0.0, '') if (x == y or (isinstance(x, str) and isinstance(y, str) and norm(x) == norm(y))) else (None, f'{path}: {str(x)[:40]!r} vs {str(y)[:40]!r}')


fc, fr = files(C), files(R)
rows = []
for rel_path in sorted(set(fc) | set(fr)):
    base = os.path.basename(rel_path)
    if SKIP_DIR.search(rel_path) or SKIP_FILE.search(base):
        rows.append((rel_path, 'SKIP', 'not compared numerically (work, logs, provenance, graphics, stub)', ''))
        continue
    if rel_path not in fr:
        rows.append((rel_path, 'MISSING', 'in control, not in run', ''))
        continue
    if rel_path not in fc:
        rows.append((rel_path, 'NEW_OK' if base in a.allow_extra else 'EXTRA', 'new output of this branch (allowed)' if base in a.allow_extra else 'in run, not in control', ''))
        continue
    p, q = fc[rel_path], fr[rel_path]
    if open(p, 'rb').read() == open(q, 'rb').read():
        rows.append((rel_path, 'IDENTICAL', '', '0'))
        continue
    ext = os.path.splitext(base)[1]
    if ext in ('.tsv', '.csv'):
        st, why, mx = cmp_table(p, q, '\t' if ext == '.tsv' else ',')
    elif ext == '.json':
        try:
            m, why = cmp_json(json.load(open(p)), json.load(open(q)))
            st, mx = ('STRUCT_DIFF' if m is None else ('IDENTICAL' if m == 0 else 'WITHIN_TOL')), (m or 0.0)
        except ValueError:
            st, why, mx = 'STRUCT_DIFF', 'not parseable as JSON', 0.0
    else:
        st, why, mx = 'STRUCT_DIFF', 'text/binary differs', 0.0
    rows.append((rel_path, st, why, f'{mx:.3g}' if mx else '0'))
with open(a.out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['file', 'status', 'detail', 'max_rel_diff'])
    w.writerows(rows)
tally = {}
for r in rows:
    tally[r[1]] = tally.get(r[1], 0) + 1
print('  '.join(f'{k}={v}' for k, v in sorted(tally.items())))
mxall = max([float(r[3]) for r in rows if r[3]] or [0.0])
print(f'largest relative difference over all compared numeric cells: {mxall:.3g} (tolerance {a.rtol:g})')
bad = [r for r in rows if r[1] in ('OUT_OF_TOL', 'STRUCT_DIFF', 'MISSING', 'EXTRA')]
for r in bad[:30]:
    print(f'  {r[1]:11s} {r[0]}  {r[2]}')
print('ACCEPTED (numeric rule)' if not bad else f'NOT ACCEPTED: {len(bad)} file(s)')
sys.exit(1 if bad else 0)
