#!/usr/bin/env python3
"""Per-sample key:value table from the GSE78220 series matrix, parsed BY KEY (characteristic rows are misaligned between samples).
usage: python3 -I gse78220_samples.py <series_matrix.txt.gz> <raw_counts.tsv.gz> <out.tsv>   (metadata only, no analysis)"""
import csv, gzip, sys

rows = {}; chars = []
for line in gzip.open(sys.argv[1], "rt"):
    if not line.startswith("!Sample_"): continue
    f = next(csv.reader([line.rstrip("\n")], delimiter="\t")); k = f[0]
    if k == "!Sample_characteristics_ch1": chars.append(f[1:])
    elif k in ("!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1", "!Sample_description"): rows[k[8:]] = f[1:]
n = len(rows["geo_accession"]); samples = []
for i in range(n):
    s = {k: v[i] for k, v in rows.items()}
    for c in chars:
        if c[i] and ":" in c[i]:
            key, val = c[i].split(":", 1)
            if key.strip() in s: sys.exit(f"duplicate key {key} in sample {i}")
            s[key.strip()] = val.strip()
    samples.append(s)
keys = list(dict.fromkeys(k for s in samples for k in s))
with open(sys.argv[3], "w") as f:
    w = csv.DictWriter(f, keys, delimiter="\t", restval="NA"); w.writeheader(); w.writerows(samples)
hdr = gzip.open(sys.argv[2], "rt").readline().rstrip("\n").split("\t")
print("samples in series matrix:", n, "| count columns:", len(hdr) - 1, "| first column:", hdr[0])
print("GSM in matrix == count columns:", sorted(s["geo_accession"] for s in samples) == sorted(hdr[1:]))
for k in ("biopsy time", "treatment", "anti-pd-1 response"):
    vals = {}
    for s in samples: vals[s.get(k, "NA")] = vals.get(s.get(k, "NA"), 0) + 1
    print(k, vals)
pt = {}
for s in samples: pt.setdefault(s.get("patient id", "NA"), []).append(s)
print("distinct patient ids:", len(pt))
for p, ss in pt.items():
    if len(ss) > 1:
        for s in ss: print("  patient", p, "|", s["geo_accession"], "|", s["title"], "|", s.get("biopsy time"), "|", s.get("anti-pd-1 response"), "|", s.get("anatomical location"))
for s in samples:
    if s.get("biopsy time") != "pre-treatment": print("  not pre-treatment:", s["geo_accession"], s["title"], s.get("patient id"), s.get("biopsy time"))
