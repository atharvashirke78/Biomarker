#!/usr/bin/env python3
"""GTF (plain or .gz) -> BED12 gene model for RSeQC (one line per transcript, exon blocks). Standard library only.
usage: gtf_to_bed12.py <annotation.gtf[.gz]> <out.bed>
Only 'exon' rows with a transcript_id are used; CDS bounds are not needed by infer_experiment.py, so thickStart/End = start/end."""
import gzip
import re
import sys


def main(gtf, out):
    tx = {}
    opener = gzip.open if gtf.endswith('.gz') else open
    with opener(gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9 or c[2] != 'exon':
                continue
            m = re.search(r'transcript_id "([^"]+)"', c[8])
            if not m:
                continue
            t = tx.setdefault(m.group(1), dict(chrom=c[0], strand=c[6], ex=[]))
            t['ex'].append((int(c[3]) - 1, int(c[4])))
    n = 0
    with open(out, 'w') as o:
        for name, t in sorted(tx.items(), key=lambda kv: (kv[1]['chrom'], min(e[0] for e in kv[1]['ex']), kv[0])):
            ex = sorted(t['ex'])
            start, end = ex[0][0], max(e[1] for e in ex)
            o.write('\t'.join(map(str, [t['chrom'], start, end, name, 0, t['strand'], start, end, 0, len(ex),
                                        ','.join(str(e[1] - e[0]) for e in ex) + ',', ','.join(str(e[0] - start) for e in ex) + ','])) + '\n')
            n += 1
    if n == 0:
        sys.exit('ERROR: no exon rows with transcript_id found in ' + gtf)
    print(f'{n} transcripts written to {out}')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
