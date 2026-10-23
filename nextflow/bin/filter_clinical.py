#!/usr/bin/env python3
"""Keep only the clinical rows whose sample is a column of the QUANT counts matrix (samples that failed QUANT QC are not in the matrix).
Same rule as `scripts/run_info.py filter-clinical` of the launcher.

usage: filter_clinical.py --clinical clinical_discover.tsv --counts counts_matrix.tsv
Writes clinical_discover.qcpass.tsv (rows kept) and excluded_after_quant.tsv (rows dropped, with the reason) into the current folder."""
import argparse

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--clinical', required=True)
ap.add_argument('--counts', required=True)
a = ap.parse_args()

cl = pd.read_csv(a.clinical, sep='\t', dtype=str, keep_default_na=False)
cols = open(a.counts).readline().rstrip('\n').split('\t')[1:]
keep = cl[cl.sample_id.isin(cols)]
drop = cl[~cl.sample_id.isin(cols)]
keep.to_csv('clinical_discover.qcpass.tsv', sep='\t', index=False)
drop.assign(reason='not in the QUANT contract counts (QC FAIL)').to_csv('excluded_after_quant.tsv', sep='\t', index=False)
print(f'clinical rows kept {len(keep)}, excluded after QUANT {len(drop)}' + (f': {", ".join(drop.sample_id)}' if len(drop) else ''))
