#!/usr/bin/env python3
"""After quantification: list the samples excluded by the sample QC (status FAIL) and stop with a clear message if too few samples are left.
usage: check_qc_min.py --qc sample_qc.tsv --min-total 3 [--min-per-group 3]
Reads sample_qc.tsv (columns srr, group, status, reason). Samples with status FAIL are excluded downstream (WARN is kept).
Always writes qc_excluded.tsv (the excluded samples and why; header only if none). Exit 2 (message on stderr) if
  - fewer than --min-total samples remain (the cohort QC needs at least 3 samples to estimate variance and a second principal component), or
  - --min-per-group > 0 and a group has fewer remaining samples than that (discover mode: analysis_min_samples_per_group, PROVISIONAL)."""
import argparse
import csv
import sys

ap = argparse.ArgumentParser()
ap.add_argument('--qc', required=True)
ap.add_argument('--min-total', type=int, default=3)
ap.add_argument('--min-per-group', type=int, default=0)
a = ap.parse_args()
rows = list(csv.DictReader(open(a.qc), delimiter='\t'))
bad = [r for r in rows if r['status'] == 'FAIL']
keep = [r for r in rows if r['status'] != 'FAIL']
with open('qc_excluded.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['sample', 'group', 'status', 'reason'])
    for r in bad:
        w.writerow([r['srr'], r['group'], r['status'], r['reason']])
excl = '; '.join(f"{r['srr']} ({r['group']}): {r['reason']}" for r in bad) or 'none'
errs = []
if len(keep) < a.min_total:
    errs.append(f'only {len(keep)} of {len(rows)} sample(s) pass the sample QC; at least {a.min_total} are needed')
if a.min_per_group:
    for g in sorted({r['group'] for r in rows}):
        n = sum(r['group'] == g for r in keep)
        if n < a.min_per_group:
            errs.append(f"group '{g}' has {n} sample(s) after QC; at least {a.min_per_group} are needed (analysis_min_samples_per_group, PROVISIONAL)")
print(f'sample QC: {len(keep)} of {len(rows)} sample(s) kept; excluded: {excl}')
if errs:
    print('RUN STOPPED after quantification (nothing downstream was run on too few samples): ' + ' | '.join(errs) + f'. Excluded by QC: {excl}', file=sys.stderr)
    sys.exit(2)
