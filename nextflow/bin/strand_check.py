#!/usr/bin/env python3
"""Second, independent strandedness estimate: RSeQC infer_experiment.py output versus the Salmon library type that set featureCounts -s.
FLAG ONLY: it never changes a count, a QC status or an exclusion.
usage: strand_check.py <sample_id> <infer_experiment.txt> <salmon_quant_dir> <stranded_min> <unstranded_max_diff>
Prints one TSV row: id, rseqc_forward, rseqc_reverse, rseqc_failed, rseqc_call, salmon_format, salmon_call, featurecounts_s, agree, note.
Cut-offs (stranded >= 0.8, unstranded when |forward - reverse| <= 0.1) are the nf-core/rnaseq values, applied to the share among determinable reads (our choice; nf-core's exact formula was not checked). Both are set in conf/thresholds.config.
forward = featureCounts -s 1 (read 1 on the transcript strand), reverse = -s 2, unstranded = -s 0."""
import json
import os
import re
import sys

COLS = ['sample', 'rseqc_forward', 'rseqc_reverse', 'rseqc_failed', 'rseqc_call', 'salmon_format', 'salmon_call', 'featurecounts_s', 'agree', 'note']


def parse_infer(text):
    """(forward, reverse, failed) fractions from infer_experiment.py output (paired or single end)."""
    def val(pat):
        m = re.search(pat + r'[^:\n]*:\s*([0-9.]+)', text)
        return float(m.group(1)) if m else None
    failed = val(r'Fraction of reads failed to determine')
    if 'PairEnd' in text:
        fwd = val(r'Fraction of reads explained by "1\+\+,1--,2\+-,2-\+"')
        rev = val(r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--"')
    else:
        fwd = val(r'Fraction of reads explained by "\+\+,--"')
        rev = val(r'Fraction of reads explained by "\+-,-\+"')
    return fwd, rev, failed


def call(fwd, rev, stranded_min, unstr_diff):
    """Cut-offs apply to the share among reads whose strand could be determined (forward + reverse = 1): typically 10-25 % of sampled reads overlap no
    unambiguous gene and are reported as 'failed to determine', which would otherwise push a clearly stranded library below the cut-off."""
    if fwd is None or rev is None or fwd + rev <= 0:
        return 'unknown'
    fwd, rev = fwd / (fwd + rev), rev / (fwd + rev)
    if fwd >= stranded_min:
        return 'forward'
    if rev >= stranded_min:
        return 'reverse'
    if abs(fwd - rev) <= unstr_diff:
        return 'unstranded'
    return 'ambiguous'


def salmon_call(fmt):
    """Same mapping as modules/featurecounts.nf."""
    if fmt in ('ISF', 'OSF', 'MSF'):
        return 'forward', 1
    if fmt in ('ISR', 'OSR', 'MSR'):
        return 'reverse', 2
    return 'unstranded', 0


def row(sid, infer_text, salmon_json_text, stranded_min=0.8, unstr_diff=0.1):
    fwd, rev, failed = parse_infer(infer_text)
    m = re.search(r'"expected_format":\s*"([A-Z]+)"', salmon_json_text)
    fmt = m.group(1) if m else 'NA'
    scall, s = salmon_call(fmt) if m else ('unknown', 'NA')
    rcall = call(fwd, rev, stranded_min, unstr_diff)
    if rcall in ('unknown', 'ambiguous') or scall == 'unknown':
        agree, note = 'undetermined', 'RSeQC call %s: not enough evidence to confirm or contradict the Salmon call (can also indicate genomic DNA contamination)' % rcall
    elif rcall == scall:
        agree, note = 'yes', ''
    else:
        agree, note = 'NO', 'RSeQC says %s but featureCounts used -s %s (%s from Salmon %s): check the library type; counts may be wrong' % (rcall, s, scall, fmt)
    f = lambda x: 'NA' if x is None else '%.4f' % x
    return [sid, f(fwd), f(rev), f(failed), rcall, fmt, scall, str(s), agree, note]


if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    sid, inf, q, smin, udiff = sys.argv[1:]
    with open(inf) as f:
        it = f.read()
    sj = ''
    p = os.path.join(q, 'lib_format_counts.json')
    if os.path.exists(p):
        with open(p) as f:
            sj = f.read()
    print('\t'.join(row(sid, it, sj, float(smin), float(udiff))))
