"""Unit tests for nextflow/bin/strand_check.py and gtf_to_bed12.py (standard library only).  python3 -m unittest tests/unit/test_strand_check.py"""
import gzip
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'nextflow', 'bin'))
import gtf_to_bed12  # noqa: E402
import strand_check as sc  # noqa: E402

PAIRED = '''This is PairEnd Data
Fraction of reads failed to determine: 0.0172
Fraction of reads explained by "1++,1--,2+-,2-+": %s
Fraction of reads explained by "1+-,1-+,2++,2--": %s
'''
SINGLE = '''This is SingleEnd Data
Fraction of reads failed to determine: 0.0100
Fraction of reads explained by "++,--": 0.0300
Fraction of reads explained by "+-,-+": 0.9600
'''
def salmon(fmt): return '{"expected_format": "%s", "compatible_fragment_ratio": 1.0}' % fmt


class TestStrand(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(sc.parse_infer(PAIRED % ('0.9000', '0.0828')), (0.9, 0.0828, 0.0172))
        self.assertEqual(sc.parse_infer(SINGLE)[:2], (0.03, 0.96))

    def test_agree_forward_and_reverse_and_unstranded(self):
        self.assertEqual(sc.row('a', PAIRED % ('0.95', '0.03'), salmon('ISF'))[8], 'yes')
        self.assertEqual(sc.row('a', PAIRED % ('0.02', '0.96'), salmon('ISR'))[8], 'yes')
        self.assertEqual(sc.row('a', PAIRED % ('0.49', '0.49'), salmon('IU'))[8], 'yes')

    def test_disagreement_is_flagged(self):
        r = sc.row('a', PAIRED % ('0.02', '0.96'), salmon('IU'))       # reverse-stranded data counted as unstranded
        self.assertEqual(r[8], 'NO')
        self.assertIn('featureCounts used -s 0', r[9])

    def test_ambiguous_and_missing_are_undetermined_not_failures(self):
        self.assertEqual(sc.row('a', PAIRED % ('0.65', '0.30'), salmon('ISF'))[8], 'undetermined')
        self.assertEqual(sc.row('a', 'garbage', salmon('ISF'))[8], 'undetermined')
        self.assertEqual(sc.row('a', PAIRED % ('0.9', '0.05'), '')[8], 'undetermined')

    def test_failed_reads_do_not_hide_a_stranded_library(self):
        r = sc.row('a', PAIRED % ('0.0063', '0.7931'), salmon('ISR'))        # 20 % undetermined reads, as in the real demo samples
        self.assertEqual((r[4], r[8]), ('reverse', 'yes'))

    def test_thresholds_come_from_arguments(self):
        self.assertEqual(sc.row('a', PAIRED % ('0.72', '0.28'), salmon('ISF'), 0.7)[4], 'forward')
        self.assertEqual(sc.row('a', PAIRED % ('0.72', '0.28'), salmon('ISF'), 0.8)[4], 'ambiguous')

    def test_gtf_to_bed12(self):
        gtf = ('chr1\tH\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
               'chr1\tH\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
               'chr1\tH\tgene\t101\t400\t.\t+\t.\tgene_id "G1";\n')
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'a.gtf.gz')
            with gzip.open(p, 'wt') as f:
                f.write(gtf)
            out = os.path.join(d, 'o.bed')
            gtf_to_bed12.main(p, out)
            c = open(out).read().strip().split('\t')
        self.assertEqual(c[:6], ['chr1', '100', '400', 'T1', '0', '+'])
        self.assertEqual(c[9:], ['2', '100,100,', '0,200,'])


if __name__ == '__main__':
    unittest.main()
