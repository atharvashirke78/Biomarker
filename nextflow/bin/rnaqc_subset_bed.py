#!/usr/bin/env python3
"""Seeded random subset of a BED12 gene model for the slow RSeQC tools (TIN, gene-body coverage). Standard library only.
usage: rnaqc_subset_bed.py <gene_model.bed12> <out.bed12> <n_transcripts> <seed> [min_exons=4]
Keeps transcripts with at least min_exons exons on nuclear chromosomes (not chrM, not alt contigs), sorted by position."""
import random
import sys


def main(src, out, n, seed, min_exons=4):
    keep = []
    with open(src) as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) >= 12 and c[0] != 'chrM' and int(c[9]) >= min_exons:
                keep.append(line)
    if not keep:
        sys.exit('ERROR: no transcript with >= %d exons in %s' % (min_exons, src))
    random.Random(seed).shuffle(keep)
    sel = sorted(keep[:n], key=lambda l: (l.split('\t')[0], int(l.split('\t')[1])))
    with open(out, 'w') as o:
        o.writelines(sel)
    print('%d of %d eligible transcripts written to %s' % (len(sel), len(keep), out))


if __name__ == '__main__':
    if len(sys.argv) not in (5, 6):
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), *(int(x) for x in sys.argv[5:]))
