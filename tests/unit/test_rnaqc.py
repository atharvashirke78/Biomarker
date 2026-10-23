"""Unit tests for nextflow/bin/rnaqc_row.py, rnaqc_flags.py, rnaqc_subset_bed.py (standard library only)."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'nextflow', 'bin'))
import rnaqc_flags as fl  # noqa: E402
import rnaqc_row as rr  # noqa: E402
import rnaqc_subset_bed as sb  # noqa: E402

RD = """Total Reads                   191228
Total Tags                    254402
Total Assigned Tags           254015
Group               Total_bases         Tag_count           Tags/Kb
CDS_Exons           211619547           250310              1.18
5'UTR_Exons         0                   0                   0.00
3'UTR_Exons         0                   0                   0.00
Introns             1865092355          3602                0.00
"""
DUP = "Occurrence\tUniqReadNumber\n1\t110753\n2\t12064\n3\t2702\n"
TIN = "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nx.bam\t30.0\t27.03\t16.0\n"
SYM = "gene_id\tgene_name\tgene_type\nG1\tRN1\trRNA\nG2\tMT-RNR1\tMt_rRNA\nG3\tA\tprotein_coding\nG4\tMT-ND1\tprotein_coding\n"
FC = "# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tx.bam\nG1\tchr1;chr1\t1\t2\t+\t10\t10\nG2\tchrM\t1\t2\t+\t10\t5\nG3\tchr2\t1\t2\t+\t10\t80\nG4\tchrM\t1\t2\t+\t10\t5\n"


class T(unittest.TestCase):
    def test_distribution(self):
        ex, intr = rr.read_distribution(RD)
        self.assertAlmostEqual(ex, 100 * 250310 / 254402, 4)
        self.assertAlmostEqual(intr, 100 * 3602 / 254402, 4)

    def test_dup_tin(self):
        tot = 110753 + 2 * 12064 + 3 * 2702
        self.assertAlmostEqual(rr.dup_rate(DUP), 100 * (1 - 110753 / tot), 4)
        self.assertEqual(rr.tin_median(TIN), 27.03)

    def test_bias(self):
        flat = 'x\n' + '\t'.join(['s'] + ['1'] * 100)
        self.assertEqual(rr.bias(flat), 1.0)
        three = 'x\n' + '\t'.join(['s'] + ['1'] * 80 + ['3'] * 20)
        self.assertEqual(rr.bias(three), 3.0)
        self.assertIsNone(rr.bias('x\ns\t1\t2'))

    def test_counts_fractions(self):
        r, m = rr.counts_fractions(FC, SYM)
        self.assertAlmostEqual(r, 15.0)       # G1 + G2 = 15 of 100
        self.assertAlmostEqual(m, 10.0)       # G2 + G4 on chrM

    def test_row_handles_missing_pieces(self):
        r = rr.row('s', '', '', '', '', FC, SYM)
        self.assertEqual(r[1:6], ['NA'] * 5)

    def test_flags(self):
        rows = [dict(sample='s%d' % i, exonic_pct='80', intronic_pct='5', dup_rate_pct='10', tin_median='30', bias_3p_5p='1.0', rrna_pct='1', mito_pct='2') for i in range(8)]
        for i, r in enumerate(rows):
            r['tin_median'] = str(30 + (i % 3))           # some spread so MAD > 0
            r['dup_rate_pct'] = str(10 + (i % 3))
        rows[0]['tin_median'] = '5'                          # degraded sample
        rows[1]['dup_rate_pct'] = '60'
        out = fl.flag(rows, 3.0, 5)
        self.assertTrue(any('tin_median low' in x for x in out[0]))
        self.assertTrue(any('dup_rate_pct high' in x for x in out[1]))
        self.assertEqual(out[4], [])
        self.assertTrue(all(v == [] for v in fl.flag(rows[:3], 3.0, 5).values()))     # too few samples: no flags

    def test_subset_is_seeded_and_filters(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, 'a.bed')
            with open(src, 'w') as f:
                for i in range(50):
                    f.write('chr1\t%d\t%d\tT%d\t0\t+\t0\t0\t0\t%d\t10,\t0,\n' % (i * 100, i * 100 + 50, i, 2 + i % 6))
                f.write('chrM\t0\t50\tM\t0\t+\t0\t0\t0\t9\t10,\t0,\n')
            a, b = os.path.join(d, 'a1'), os.path.join(d, 'a2')
            sb.main(src, a, 10, 7)
            sb.main(src, b, 10, 7)
            self.assertEqual(open(a).read(), open(b).read())
            self.assertNotIn('chrM', open(a).read())
            self.assertTrue(all(int(l.split('\t')[9]) >= 4 for l in open(a)))


if __name__ == '__main__':
    unittest.main()
