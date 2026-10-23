#!/usr/bin/env python3
"""One per-sample row of RNA-QC metrics from RSeQC outputs and the featureCounts file. FLAG ONLY: never changes a count or a QC status.
usage: rnaqc_row.py <sample_id> <read_distribution.txt> <dup.pos.DupRate.xls> <tin.summary.txt> <geneBodyCoverage.txt|NONE> <featureCounts.txt> <gene_id_to_symbol.tsv>
Columns: sample, exonic_pct, intronic_pct, dup_rate_pct, tin_median, bias_3p_5p, rrna_pct, mito_pct, note
- exonic_pct / intronic_pct: share of RSeQC 'Total Tags' in exons (CDS+UTR groups) and introns (a rough proxy for gDNA or pre-mRNA).
- dup_rate_pct: position-based duplication, 100 * (1 - reads at positions seen once / all reads) from read_duplication.py.
- tin_median: median transcript integrity number on a seeded random subset of transcripts (lower = more degradation).
- bias_3p_5p: mean gene-body coverage in percentile bins 81-100 divided by bins 1-20 (own definition; transcripts are oriented 5'->3'; >1 = 3' bias). Not Picard's metric.
- rrna_pct, mito_pct: share of featureCounts-assigned reads on genes of type rRNA / Mt_rRNA, and on chrM genes."""
import csv
import re
import sys

COLS = ['sample', 'exonic_pct', 'intronic_pct', 'dup_rate_pct', 'tin_median', 'bias_3p_5p', 'rrna_pct', 'mito_pct', 'note']


def read_distribution(text):
    tags = re.search(r'Total Tags\s+(\d+)', text)
    if not tags or int(tags.group(1)) == 0:
        return None, None
    total = int(tags.group(1))
    g = {}
    for line in text.splitlines():
        c = line.split()
        if len(c) == 4 and c[0] in ("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns"):
            g[c[0]] = int(c[2])
    ex = g.get('CDS_Exons', 0) + g.get("5'UTR_Exons", 0) + g.get("3'UTR_Exons", 0)
    return 100.0 * ex / total, 100.0 * g.get('Introns', 0) / total


def dup_rate(text):
    occ = [l.split() for l in text.splitlines()[1:] if len(l.split()) == 2]
    total = sum(int(a) * int(b) for a, b in occ)
    if total == 0:
        return None
    once = sum(int(b) for a, b in occ if a == '1')
    return 100.0 * (1 - once / total)


def tin_median(text):
    for line in text.splitlines()[1:]:
        c = line.split('\t')
        if len(c) >= 3:
            try:
                return float(c[2])
            except ValueError:
                return None
    return None


def bias(text):
    if not text:
        return None
    rows = [l.split('\t') for l in text.strip().splitlines() if l.strip()]
    if len(rows) < 2:
        return None
    v = [float(x) for x in rows[1][1:]]
    if len(v) != 100:
        return None
    five, three = sum(v[0:20]) / 20, sum(v[80:100]) / 20
    return three / five if five > 0 else None


def counts_fractions(fc_text, symbols_text):
    rrna_ids = {r['gene_id'] for r in csv.DictReader(symbols_text.splitlines(), delimiter='\t') if r.get('gene_type') in ('rRNA', 'Mt_rRNA')}
    total = rrna = mito = 0
    lines = [l for l in fc_text.splitlines() if not l.startswith('#')]
    for line in lines[1:]:
        c = line.split('\t')
        if len(c) < 7:
            continue
        n = int(c[6])
        total += n
        if c[0] in rrna_ids:
            rrna += n
        if c[1].split(';')[0] == 'chrM':
            mito += n
    if total == 0:
        return None, None
    return 100.0 * rrna / total, 100.0 * mito / total


def row(sid, rd, dup, tin, gbc, fc, sym):
    f = lambda x, d=2: 'NA' if x is None else ('%.' + str(d) + 'f') % x
    ex, intr = read_distribution(rd)
    rr, mt = counts_fractions(fc, sym)
    notes = []
    b = bias(gbc)
    if gbc and b is None:
        notes.append('gene-body coverage unreadable')
    return [sid, f(ex), f(intr), f(dup_rate(dup)), f(tin_median(tin)), f(b, 3), f(rr, 3), f(mt, 3), '; '.join(notes)]


if __name__ == '__main__':
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    sid, a, b_, c, d, e, g = sys.argv[1:]
    rd = lambda p: open(p).read() if p != 'NONE' else ''
    print('\t'.join(row(sid, rd(a), rd(b_), rd(c), rd(d), rd(e), rd(g))))
