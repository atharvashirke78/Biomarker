#!/usr/bin/env python3
# usage: merge_featurecounts.py <dir with <label>.featureCounts.txt> <out.tsv> <sample_map.tsv>   (integer genes x samples matrix; columns = SRR)
import csv, glob, os, sys
d, out, smap = sys.argv[1:4]
lab2srr = {}
for i, l in enumerate(open(smap)):
    if i: p = l.rstrip('\n').split('\t'); lab2srr[p[1]] = p[0]
cols, genes = {}, None
for f in sorted(glob.glob(os.path.join(d, "*.featureCounts.txt"))):
    lab = os.path.basename(f).replace(".featureCounts.txt", "")
    if lab not in lab2srr: sys.exit("ERROR: no SRR for label " + lab)
    s = lab2srr[lab]
    ids, vals = [], []
    with open(f) as fh:
        for l in fh:
            if l.startswith("#") or l.startswith("Geneid"): continue
            p = l.rstrip("\n").split("\t"); ids.append(p[0]); vals.append(int(p[-1]))
    if genes is None: genes = ids
    elif ids != genes: sys.exit("ERROR: gene order differs in " + f)
    cols[s] = vals
names = sorted(cols)
with open(out, "w") as o:
    w = csv.writer(o, delimiter="\t"); w.writerow(["gene_id"] + names)
    for i, g in enumerate(genes): w.writerow([g] + [cols[n][i] for n in names])
print("merged %d genes x %d samples" % (len(genes), len(names)))
