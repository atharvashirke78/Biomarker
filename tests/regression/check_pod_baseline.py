#!/usr/bin/env python3
"""Two-level regression, level 2: byte comparison of a regression run with the pod baseline (read-only).

usage (project root): env/prefix/pyqc/bin/python3 -I tests/regression/check_pod_baseline.py --baseline <run dir> --run <run dir> --out <new tsv>
Both dirs come from tests/regression/run_regression.sh. Every file under the baseline's published output folders is compared:
  DATA  (.tsv .csv .json from analysis steps, .txt model/spec files, .png, .rds, everything not listed below): byte-identical required.
  META  (*.provenance.txt, RESULTS_INDEX*.md, DISCOVER_SUMMARY.md, MODEL_MANIFEST.tsv, manifest_*.json, run_info*, tool_versions.txt,
         sessionInfo.txt): compared after masking the run folder path, Nextflow run names, timestamps and durations only
         (file hashes are NOT masked). Masked tokens are counted.
  SKIP  (Nextflow work/, .nextflow*, logs, traces, html reports, timeline): listed, not compared.
Exit 0 only if every DATA and META file is IDENTICAL and no file is missing or extra. Nothing is changed to make a comparison match."""
import argparse
import os
import re
import sys

ap = argparse.ArgumentParser()
ap.add_argument('--baseline', required=True)
ap.add_argument('--run', required=True)
ap.add_argument('--out', required=True)
a = ap.parse_args()
if os.path.exists(a.out):
    sys.exit('ERROR: output exists (never overwritten): ' + a.out)
B, N = os.path.realpath(a.baseline), os.path.realpath(a.run)

SKIP_DIRS = {'work', '.nextflow'}
SKIP_FILE = re.compile(r'(\.log(\.\d+)?$|^\.nextflow|trace.*\.tsv$|report.*\.html$|timeline.*\.html$|^exit_code\.txt$|\.console$|^run_info\.txt$|^diagnostic_exit_code\.txt$|^diagnostic\.log$)')
META_FILE = re.compile(r'(\.provenance\.txt$|^RESULTS_INDEX.*\.md$|^DISCOVER_SUMMARY\.md$|^MODEL_MANIFEST\.tsv$|^manifest_.*\.json$|^run_info.*$|^tool_versions\.txt$|^sessionInfo\.txt$|^command\.txt$|^regression_summary\.tsv$|^rule4_summary\.json$)')
TIME = [re.compile(p) for p in (r'\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(\.\d+)?Z?', r'\b\d+(\.\d+)?\s?(ms|s|m|h)\b(?=[\s,;)]|$)',
                                r'\b(Mon|Tue|Wed|Thu|Fri|Sat|Sun) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) [ \d]\d \d{2}:\d{2}:\d{2} \w+ \d{4}')]


def run_names(root):
    names = set()
    for dp, dn, fs in os.walk(root):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        for f in fs:
            if f.startswith('RESULTS_INDEX') or f == 'DISCOVER_SUMMARY.md':
                head = open(os.path.join(dp, f), errors='replace').readline()
                names.update(re.findall(r'(?:run: |summary: )([a-z]+_[a-z]+)', head))
    return names


def mask(text, root, names):
    n = 0
    text, k = re.subn(re.escape(root), '<RUN>', text); n += k
    text, k = re.subn(r'/work/[0-9a-f]{2}/[0-9a-f]{30}/', '/work/<TASK>/', text); n += k   # Nextflow task folder hash
    # the ONLY masked file hash: DISCOVER_SUMMARY.md contains the random Nextflow run name, so its sha256 in PACKAGE_MODEL.provenance.txt
    # cannot match between runs (the summary itself is compared as META above, with the run name masked)
    text, k = re.subn(r'(?m)^[0-9a-f]{64}(?=  DISCOVER_SUMMARY\.md$)', '<HASH_OF_RUN_NAMED_FILE>', text); n += k
    for nm in names:
        text, k = re.subn(r'\b' + re.escape(nm) + r'\b', '<RUNNAME>', text); n += k
    for rx in TIME:
        text, k = rx.subn('<TIME>', text); n += k
    return text, n


nb, nn = run_names(B), run_names(N)
rows = []


def walk(root):
    out = {}
    for dp, dn, fs in os.walk(root):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        for f in fs:
            p = os.path.join(dp, f)
            out[os.path.relpath(p, root)] = p
    return out


fb, fn = walk(B), walk(N)
for rel in sorted(set(fb) | set(fn)):
    base = os.path.basename(rel)
    if SKIP_FILE.search(base):
        rows.append((rel, 'SKIP', 'not compared (log/trace/report)'))
        continue
    if rel not in fn:
        rows.append((rel, 'MISSING', 'in baseline, not in run'))
        continue
    if rel not in fb:
        rows.append((rel, 'EXTRA', 'in run, not in baseline'))
        continue
    pb, pn = fb[rel], fn[rel]
    if os.path.islink(pb) or os.path.islink(pn):
        ok = os.path.islink(pb) == os.path.islink(pn) and os.path.basename(os.readlink(pb) if os.path.islink(pb) else '') == os.path.basename(os.readlink(pn) if os.path.islink(pn) else '')
        rows.append((rel, 'IDENTICAL' if ok else 'DIFFERENT', 'symlink'))
        continue
    xb, xn = open(pb, 'rb').read(), open(pn, 'rb').read()
    if xb == xn:
        rows.append((rel, 'IDENTICAL', 'bytes'))
        continue
    if META_FILE.search(base):
        tb, kb = mask(xb.decode('utf-8', 'replace'), B, nb)
        tn, kn = mask(xn.decode('utf-8', 'replace'), N, nn)
        if tb == tn:
            rows.append((rel, 'IDENTICAL', f'META after masking ({kb} / {kn} tokens masked)'))
            continue
        lb, ln = tb.splitlines(), tn.splitlines()
        i = next((i for i in range(min(len(lb), len(ln))) if lb[i] != ln[i]), min(len(lb), len(ln)))
        rows.append((rel, 'DIFFERENT', f'META; first differing line {i + 1}: base={lb[i][:160] if i < len(lb) else "<eof>"!r} run={ln[i][:160] if i < len(ln) else "<eof>"!r}'))
        continue
    try:
        lb, ln = xb.decode().splitlines(), xn.decode().splitlines()
        i = next((i for i in range(min(len(lb), len(ln))) if lb[i] != ln[i]), min(len(lb), len(ln)))
        nd = sum(1 for x, y in zip(lb, ln) if x != y) + abs(len(lb) - len(ln))
        rows.append((rel, 'DIFFERENT', f'DATA; {nd} lines differ; first {i + 1}: base={lb[i][:160] if i < len(lb) else "<eof>"!r} run={ln[i][:160] if i < len(ln) else "<eof>"!r}'))
    except UnicodeDecodeError:
        rows.append((rel, 'DIFFERENT', f'DATA binary; sizes {len(xb)} vs {len(xn)}'))

with open(a.out, 'w') as fh:
    fh.write('file\tstatus\tdetail\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
from collections import Counter
c = Counter(r[1] for r in rows)
print(f'baseline {B}\nrun      {N}\n' + ', '.join(f'{k} {v}' for k, v in sorted(c.items())))
for r in rows:
    if r[1] in ('DIFFERENT', 'MISSING', 'EXTRA'):
        print('\t'.join(r)[:400])
sys.exit(0 if not any(r[1] in ('DIFFERENT', 'MISSING', 'EXTRA') for r in rows) else 1)
