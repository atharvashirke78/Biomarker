#!/usr/bin/env python3
"""Independent check of both workbooks: csv-module recomputation (no shared code with build_comparison.py), formulas vs values copy, category sums, causes."""
import csv, re, sys
import pandas as pd, openpyxl
B = "/workspace/Biomarker_pipeline"; R = B + "/results"; V = R + "/09_verification"; X = V + "/input/MOESM2_sheets"; D = V + "/step11_comparison_workbook/"
F_LIVE = D + "Hosseini2023_vs_Pipeline_comparison.xlsx"; F_VAL = D + "Hosseini2023_vs_Pipeline_comparison_values.xlsx"
want = ["Summary", "Gene_overlap", "Missed_breakdown", "Opposite_genes", "Not_in_table_causes", "Pathways_GO_KEGG", "Readcounts_S2", "Gene18", "Method_notes"]
fail = []
def chk(ok, msg):
    print(("OK   " if ok else "FAIL ") + msg); (None if ok else fail.append(msg))
# a) reopen both
wl = openpyxl.load_workbook(F_LIVE); wv = openpyxl.load_workbook(F_VAL)
for nm, w in (("live", wl), ("values", wv)):
    chk(w.sheetnames == want, f"{nm}: sheets {w.sheetnames}")
    errs = [(ws.title, c.coordinate) for ws in w for row in ws.iter_rows() for c in row if isinstance(c.value, str) and re.match(r"^#(REF|DIV/0|VALUE|NAME|N/A|NUM|NULL)", c.value)]
    nf = sum(1 for ws in w for row in ws.iter_rows() for c in row if isinstance(c.value, str) and c.value.startswith("="))
    chk(not errs, f"{nm}: no error values ({errs}); formula cells: {nf}")
    if nm == "values": chk(nf == 0, "values copy has no formulas")
    else: chk(nf > 0, "live file has formulas")
for f in (F_LIVE, F_VAL): sh = pd.read_excel(f, sheet_name=None, engine="openpyxl"); print("  pandas read:", f.split("/")[-1], {k: v.shape for k, v in sh.items()})
# b) independent recomputation
def de(path):
    best = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        g = r["gene_name"]
        if not g: continue
        padj = None if r["padj"] in ("NA", "") else float(r["padj"]); key = (padj is not None, float(r["baseMean"]))
        if g not in best or key > best[g][0]: best[g] = (key, float(r["log2FoldChange"]), padj)
    return {g: (l, p, p is not None and p < 0.05 and abs(l) >= 1) for g, (k, l, p) in best.items()}
paper = {r["gene"].strip(): r["direction_in_cancer"] for r in csv.DictReader(open(V + "/input/hosseini_2023_TableS3_genes.tsv"), delimiter="\t")}
def counts(t):
    o = {}
    for d in ("up", "down"):
        pl = {g for g, x in paper.items() if x == d}; ours = {g for g, v in t.items() if v[2] and (v[0] > 0) == (d == "up")}
        m = {g for g in pl if g in t and t[g][2] and (t[g][0] > 0) == (d == "up")}; opp = {g for g in pl if g in t and t[g][2]} - m
        o[d] = (len(pl), len(ours), len(m), len(opp), len(pl) - len(m) - len(opp), len(ours - set(paper)))
    return o
fcT = de(R + "/06_dge/deseq2_results_all.tsv"); saT = de(R + "/06_dge/deseq2_salmon_results_all.tsv"); fcc = counts(fcT); sac = counts(saT)
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
exp = {"Down genes (featureCounts, primary)": fcc["down"], "Up genes (featureCounts, primary)": fcc["up"], "All DEGs (featureCounts, primary)": add(fcc["up"], fcc["down"]),
       "Down genes (Salmon, secondary)": sac["down"], "Up genes (Salmon, secondary)": sac["up"], "All DEGs (Salmon, secondary)": add(sac["up"], sac["down"])}
def norm(x): return re.sub(r"[^a-z0-9]+", " ", x.lower()).strip()
def paper_terms(f, kegg):
    out = set()
    for r in list(csv.reader(open(f), delimiter="\t", quoting=csv.QUOTE_NONE))[2:]:
        if len(r) > 4 and float(r[4]) < 0.05: out.add(norm(r[1]) if kegg else re.search(r"\((GO:\d+)\)", r[1]).group(1))
    return out
def our_terms(f, kegg): return {(norm(r["Description"]) if kegg else r["ID"]) for r in csv.DictReader(open(f), delimiter="\t") if float(r["p.adjust"]) < 0.05}
ora = R + "/07_pathway/ora/"
def pc(typ):
    kegg = typ == "KEGG"; P = {}; O = {}
    for d in ("up", "down"):
        P[d] = paper_terms(X + ("/Table_S6_Upregulated_.tsv" if d == "up" else "/Table_S6_Downregulated_.tsv") if kegg else X + ("/Table_S5_Biological_process_.tsv" if d == "up" else "/Table_S4_Biological_process_.tsv"), kegg)
        O[d] = our_terms(ora + f"ora_{d}_{'KEGG' if kegg else 'GO_BP'}.tsv", kegg)
    res = {}
    for d in ("up", "down"):
        e = "down" if d == "up" else "up"; m = P[d] & O[d]; opp = (P[d] & O[e]) - O[d]; ex = O[d] - P[d] - P[e]
        res[d] = (len(P[d]), len(O[d]), len(m), len(opp), len(P[d]) - len(m) - len(opp), len(ex))
    return res
go = pc("GO"); ke = pc("KEGG")
exp["GO BP terms (up + down, padj < 0.05)"] = add(go["up"], go["down"]); exp["KEGG down pathways (padj < 0.05)"] = ke["down"]; exp["KEGG up pathways (padj < 0.05)"] = ke["up"]
for nm, w in (("live", wl), ("values", wv)):
    ws = w["Summary"]; seen = set()
    for row in ws.iter_rows(min_row=2, max_row=14, values_only=True):
        if row[0] in exp:
            got = tuple(row[1:7]); chk(got == tuple(exp[row[0]]), f"{nm} Summary {row[0]}: {got} vs independent {exp[row[0]]}"); seen.add(row[0])
            chk(got[2] + got[3] + got[4] == got[0], f"   adds up to paper total ({got[2]}+{got[3]}+{got[4]} = {got[0]})")
    chk(seen == set(exp), f"{nm} Summary: all 9 gene/term rows found")
    r18 = [r for r in ws.iter_rows(min_row=2, max_row=14, values_only=True) if str(r[0]).startswith("18-gene")][0]
    print("  18-gene row:", r18[:7]); chk(r18[2] == "significant 12, expected direction 17/18" and r18[3] == 12 and r18[2 + 2] == 0 and r18[5] == 6, f"{nm} 18-gene row text and numbers")
    chk(ws["I1"].value == "Our hits also in paper list %" and ws["E1"].value == "Opposite direction", f"{nm} Summary headers")
# c) formulas vs values copy (evaluate the live formulas independently)
sl, sv = wl["Summary"], wv["Summary"]
for r in range(2, 14):
    D_, B_, C_, E_, F_ = (sl.cell(r, c).value for c in (4, 2, 3, 5, 6))
    if sl.cell(r, 8).value and str(sl.cell(r, 8).value).startswith("="):
        assert sl.cell(r, 8).value == f'=IFERROR(D{r}/B{r},"n/a")' and sl.cell(r, 9).value == f'=IFERROR(D{r}/C{r},"n/a")' and sl.cell(r, 10).value == f'=IF(D{r}+E{r}+F{r}=B{r},"OK","CHECK")', r
        e8 = D_ / B_; e9 = (D_ / C_) if isinstance(C_, (int, float)) and C_ else "n/a"
        ok = abs(sv.cell(r, 8).value - e8) < 1e-12 and ((sv.cell(r, 9).value == e9) if e9 == "n/a" else abs(sv.cell(r, 9).value - e9) < 1e-12) and sv.cell(r, 10).value == "OK"
        chk(ok, f"values copy row {r} ({sl.cell(r, 1).value[:32]}): recall {sv.cell(r, 8).value:.4f}, hits-in-paper {sv.cell(r, 9).value if isinstance(sv.cell(r, 9).value, str) else round(sv.cell(r, 9).value, 4)}")
for nm in ("Missed_breakdown",):
    a, b = wl[nm], wv[nm]
    for r in range(2, 9):
        for c, cc in ((3, 2), (5, 4)):
            chk(abs(b.cell(r, c).value - b.cell(r, cc).value / b.cell(8, cc).value) < 1e-12 or r == 8, f"values copy Missed_breakdown ({r},{c})") if False else None
    tot = [b.cell(8, 2).value, b.cell(8, 4).value]; chk(tot == [285, 4752], f"Missed_breakdown totals (values copy) {tot}")
    chk(all(abs(b.cell(r, 3).value - b.cell(r, 2).value / 285) < 1e-12 and abs(b.cell(r, 5).value - b.cell(r, 4).value / 4752) < 1e-12 for r in range(2, 8)), "Missed_breakdown % = n / paper list (values copy)")
# d) category sums + opposite genes + causes
go_df = pd.read_excel(F_VAL, sheet_name="Gene_overlap", usecols="A:G")
pu = {g for g, x in paper.items() if x == "up"}; pdn = {g for g, x in paper.items() if x == "down"}
union = set(paper) | {g for g, v in fcT.items() if v[2]}
chk(len(go_df) == len(union) == go_df.gene.nunique(), f"Gene_overlap rows {len(go_df)} = union {len(union)}")
for d, pl in (("up", pu), ("down", pdn)):
    sub = go_df[go_df["paper direction"] == d]; chk(len(sub) == len(pl), f"paper {d}: Matched+Opposite+Paper-only rows {len(sub)} = list size {len(pl)}  {dict(sub.category.value_counts())}")
opp_exp = {g for g in paper if g in fcT and fcT[g][2] and ((fcT[g][0] > 0) != (paper[g] == "up"))}
op = pd.read_excel(F_VAL, sheet_name="Opposite_genes", usecols="A:E").dropna(subset=["gene"]); op = op[op.gene.isin(paper)]
chk(set(op.gene) == opp_exp and len(opp_exp) == 10, f"Opposite_genes sheet = independent list ({len(opp_exp)}): {sorted(opp_exp)}")
# causes: independent from annotation + count matrix
ann = {};
for r in csv.DictReader(open(B + "/ref/gene_id_to_symbol.tsv"), delimiter="\t"): ann.setdefault(r["gene_name"], []).append((r["gene_id"], r["gene_type"]))
cm = list(csv.reader(open(R + "/04_counts/counts_featurecounts.tsv"), delimiter="\t")); hdr = cm[0][1:]
used = [l for l in open(R + "/06_dge/deseq2_normalized_counts.tsv").readline().rstrip("\n").split("\t")[1:]]; idx = [hdr.index(u) for u in used]
stat = {}
for r in cm[1:]:
    v = [float(r[1 + i]) for i in idx]; stat[r[0]] = (sum(x >= 10 for x in v) >= 5, sum(v))
tab = set(fcT)
caus = {}
for g, d in paper.items():
    if g in tab: continue
    if g not in ann: c = "nosym"
    else:
        ids = [i for i, _ in ann[g]]
        if any(stat[i][0] for i in ids): c = "other"
        else:
            best = max(ids, key=lambda i: stat[i][1]); bt = dict(ann[g])[best]; c = "pc" if bt == "protein_coding" else ("lnc" if bt == "lncRNA" else "otherbio")
    caus[g] = (d, c)
ws = wv["Not_in_table_causes"]; sheet_tab = {ws.cell(r, 9).value: (ws.cell(r, 10).value, ws.cell(r, 11).value) for r in range(2, 7)}
names = {"pc": "Removed by our low-count filter: protein-coding gene", "lnc": "Removed by our low-count filter: lncRNA", "otherbio": "Removed by our low-count filter: other biotype (pseudogene, miRNA, etc.)",
         "nosym": "No matching gene symbol in our annotation (alias / renamed symbol)", "other": "Other (in annotation and passes the filter, but not in our table)"}
for k, nme in names.items():
    e = (sum(1 for d, c in caus.values() if c == k and d == "down"), sum(1 for d, c in caus.values() if c == k and d == "up")); chk(sheet_tab[nme] == e, f"cause '{nme[:55]}': sheet down/up {sheet_tab[nme]} vs independent {e}")
chk(len(caus) == 1438 and sum(1 for d, c in caus.values() if d == "up") == 1431, f"not-in-table genes: {len(caus)} (up {sum(1 for d, c in caus.values() if d == 'up')}, down {sum(1 for d, c in caus.values() if d == 'down')})")
lst = pd.read_excel(F_VAL, sheet_name="Not_in_table_causes", usecols="A:C").dropna(subset=["gene"]); lst = lst[lst.gene.isin(paper)]
rev = {v: k for k, v in names.items()}; chk(all(caus[g][1] == rev[c] for g, c in zip(lst.gene, lst.cause)) and len(lst) == 1438, "per-gene causes match the independent calculation")
rc = wv["Readcounts_S2"]; chk(sum(1 for r in rc.iter_rows(min_row=2, values_only=True) if r[3] is True) == 20, "Readcounts: 20/20 TRUE")
print("RESULT:", "ALL OK" if not fail else f"PROBLEMS: {fail}"); sys.exit(1 if fail else 0)
