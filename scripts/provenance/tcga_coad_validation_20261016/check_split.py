#!/usr/bin/env python3
"""PREREG (b)/D1: independent check that no patient is on both sides. usage: python3 -I check_split.py <clinical.tsv> <holdout_samples.txt> [<training_samples.txt>]"""
import csv, sys

cl = {r["sample_id"]: r["patient_id"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
ho = {l.strip() for l in open(sys.argv[2]) if l.strip()}
tr = {l.strip() for l in open(sys.argv[3]) if l.strip()} if len(sys.argv) > 3 else set(cl) - ho
unknown = (ho | tr) - set(cl); overlap_s = ho & tr
both = sorted({cl[s] for s in ho} & {cl[s] for s in tr})
print(f"training samples {len(tr)}, hold-out samples {len(ho)}, unknown IDs {len(unknown)}, samples on both sides {len(overlap_s)}, patients on both sides {len(both)}")
print("RESULT", "PASS" if not (unknown or overlap_s or both) and tr | ho == set(cl) else "FAIL")
