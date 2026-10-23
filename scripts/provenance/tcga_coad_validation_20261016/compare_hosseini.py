#!/usr/bin/env python3
"""PREREG (f) D2, D3, D6: our training-fit DESeq2 table vs Hosseini 2023 S3 / S7 (and S4/S5 GO BP vs our ORA). Definitions exactly as in PREREG.md.
usage: python3 -I compare_hosseini.py <run1 01_de/deseq2_results_all.tsv> <results/09_verification/input dir> <pathways dir> <out dir (must not exist)>"""
import json, os, re, sys
import pandas as pd

de_f, hin, pw, out = sys.argv[1:5]
os.makedirs(out)                                        # fails if it exists
P_CUT, LFC_CUT = 0.05, 1.0                              # rule 1 (config/candidate_rules.yaml rule1_significance)
d = pd.read_csv(de_f, sep="\t")
u = d[d.pvalue.notna() & d.gene_name.notna() & (d.gene_name.astype(str) != "")]
u = u.sort_values("baseMean", ascending=False, kind="mergesort").drop_duplicates("gene_name")     # one row per symbol: highest baseMean
u = u.assign(sig=(u.padj < P_CUT) & (u.log2FoldChange_shrunk.abs() >= LFC_CUT), ours=u.log2FoldChange_shrunk.map(lambda x: "up" if x > 0 else "down"))
s3 = pd.read_csv(os.path.join(hin, "hosseini_2023_TableS3_genes.tsv"), sep="\t"); s3["gene"] = s3.gene.str.strip(); s3 = s3.drop_duplicates("gene")
m = u.merge(s3[["gene", "direction_in_cancer"]], left_on="gene_name", right_on="gene")
sig_in = m[m.sig]; agree = int((sig_in.ours == sig_in.direction_in_cancer).sum())
D2 = agree / len(sig_in) if len(sig_in) else float("nan")
# D3
s7 = open(os.path.join(hin, "MOESM2_sheets", "Table_S7_p_value_.tsv")).read().split("\n")[1].split("\t")
s7 = [g.strip() for g in s7 if g.strip()]; assert len(s7) == 40, len(s7)
s7d = dict(zip(s3.gene, s3.direction_in_cancer)); assert all(g in s7d for g in s7)
rows = []
for g in s7:
    z = u[u.gene_name == g]
    if len(z) == 0 or pd.isna(z.padj.iloc[0]): rows.append(dict(gene=g, published=s7d[g], testable=False)); continue
    z = z.iloc[0]
    rows.append(dict(gene=g, published=s7d[g], testable=True, baseMean=z.baseMean, log2FC_shrunk=z.log2FoldChange_shrunk, padj=z.padj, rule1=bool(z.sig), ours=z.ours,
                     pass_=bool(z.sig and z.ours == s7d[g])))
t7 = pd.DataFrame(rows); t7.to_csv(os.path.join(out, "D3_S7_genes.tsv"), sep="\t", index=False)
nt = int(t7.testable.sum()); npass = int(t7.get("pass_", pd.Series(dtype=bool)).fillna(False).sum()); D3 = npass / nt if nt else float("nan")
# D6 recall / overlap
s3_in_u = m; matched = int((m.sig & (m.ours == m.direction_in_cancer)).sum()); opposite = int((m.sig & (m.ours != m.direction_in_cancer)).sum())
missed_tested = int((~m.sig).sum()); not_in_universe = int(len(s3) - len(m)); extra = int(u.sig.sum() - m.sig.sum())
by_dir = {dd: dict(s3_total=int((s3.direction_in_cancer == dd).sum()), in_universe=int((m.direction_in_cancer == dd).sum()),
                   matched=int((m.sig & (m.ours == dd) & (m.direction_in_cancer == dd)).sum())) for dd in ("up", "down")}
# D6 pathways: GO BP IDs, ours (ora.R, p.adjust < 0.05) vs S5 (up) / S4 (down), Adjusted.P.value < 0.05
def go_ids(series): return {x for s in series for x in re.findall(r"GO:\d{7}", str(s))}
pwr = {}
for dd, sheet in (("up", "Table_S5_Biological_process_.tsv"), ("down", "Table_S4_Biological_process_.tsv")):
    h = pd.read_csv(os.path.join(hin, "MOESM2_sheets", sheet), sep="\t", skiprows=1); H = go_ids(h.Term[h["Adjusted.P.value"] < 0.05])
    f = os.path.join(pw, f"ora_{dd}_GO_BP.tsv")
    o = pd.read_csv(f, sep="\t") if os.path.getsize(f) > 0 else pd.DataFrame(columns=["ID", "p.adjust"])
    O = set(o.ID[o["p.adjust"] < 0.05]) if len(o) else set()
    pwr[dd] = dict(hosseini_terms=len(H), our_terms=len(O), overlap=len(H & O), jaccard=round(len(H & O) / len(H | O), 4) if H | O else None,
                   share_of_hosseini_terms_found=round(len(H & O) / len(H), 4) if H else None)
R = {"universe_symbols": int(len(u)), "our_significant_symbols_rule1": int(u.sig.sum()),
     "D2": dict(significant_and_in_S3=int(len(sig_in)), same_direction=agree, agreement=D2, threshold=0.90, status="PASS" if D2 >= 0.90 else "FAIL"),
     "D3": dict(S7_genes=40, testable=nt, untestable=[r["gene"] for r in rows if not r["testable"]], significant_in_published_direction=npass, fraction=D3, threshold=0.75,
                status="PASS" if D3 >= 0.75 else "FAIL"),
     "D6_reported": dict(S3_genes=int(len(s3)), S3_in_universe=int(len(m)), S3_not_in_universe=not_in_universe, matched=matched, opposite=opposite,
                         S3_tested_not_significant=missed_tested, extra_ours_not_in_S3=extra, recall_vs_S3_in_universe=round(matched / len(m), 4),
                         recall_vs_all_S3=round(matched / len(s3), 4), by_direction=by_dir, pathways_GO_BP=pwr)}
json.dump(R, open(os.path.join(out, "hosseini_comparison.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
