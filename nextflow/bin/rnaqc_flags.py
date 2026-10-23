#!/usr/bin/env python3
"""Merge per-sample RNA-QC rows and add cohort-derived flags. FLAG ONLY: no sample is excluded and no QC status changes.
usage: rnaqc_flags.py <rows.tsv (no header)> <sample_map.tsv> <mad_cutoff> <min_samples> > rna_qc_metrics.tsv
A metric is flagged for a sample when its robust z-score, (x - median) / (1.4826 * MAD), is beyond +-mad_cutoff in the unfavourable direction
(high intronic %, duplication %, rRNA %, mito %, 3'/5' ratio; low exonic %, TIN). There is no published pass/fail cut-off for these metrics
(docs/review/03_qc_reproducibility.md): the flags are relative to the cohort and are OUR choice. With fewer than min_samples samples or MAD = 0 nothing is flagged."""
import statistics
import sys

METRICS = [('exonic_pct', -1), ('intronic_pct', 1), ('dup_rate_pct', 1), ('tin_median', -1), ('bias_3p_5p', 1), ('rrna_pct', 1), ('mito_pct', 1)]
HEAD = ['srr', 'label', 'exonic_pct', 'intronic_pct', 'dup_rate_pct', 'tin_median', 'bias_3p_5p', 'rrna_pct', 'mito_pct', 'flags']


def flag(rows, cutoff, min_n):
    """rows: list of dicts with metric strings. Returns {index: [flag text]}"""
    out = {i: [] for i in range(len(rows))}
    for m, direction in METRICS:
        vals = [(i, float(r[m])) for i, r in enumerate(rows) if r.get(m) not in (None, '', 'NA')]
        if len(vals) < min_n:
            continue
        med = statistics.median(v for _, v in vals)
        mad = statistics.median(abs(v - med) for _, v in vals)
        if mad == 0:
            continue
        for i, v in vals:
            z = (v - med) / (1.4826 * mad)
            if direction * z > cutoff:
                out[i].append('%s %s (robust z %.1f)' % (m, 'high' if direction > 0 else 'low', z))
    return out


def main(rows_path, map_path, cutoff, min_n):
    cols = ['sample', 'exonic_pct', 'intronic_pct', 'dup_rate_pct', 'tin_median', 'bias_3p_5p', 'rrna_pct', 'mito_pct', 'note']
    rows = [dict(zip(cols, l.rstrip('\n').split('\t'))) for l in open(rows_path) if l.strip()]
    srr = {}
    for l in list(open(map_path))[1:]:
        c = l.rstrip('\n').split('\t')
        if len(c) >= 2:
            srr[c[1]] = c[0]
    rows.sort(key=lambda r: srr.get(r['sample'], r['sample']))
    fl = flag(rows, cutoff, min_n)
    print('\t'.join(HEAD))
    for i, r in enumerate(rows):
        note = '; '.join(fl[i]) if fl[i] else ('not enough samples for cohort flags' if len(rows) < min_n else '')
        print('\t'.join([srr.get(r['sample'], r['sample']), r['sample']] + [r[m] for m, _ in METRICS] + [note]))


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]))
