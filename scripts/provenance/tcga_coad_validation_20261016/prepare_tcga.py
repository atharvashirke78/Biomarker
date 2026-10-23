#!/usr/bin/env python3
"""PREREG (a)-(c): one tumour (+ at most one normal) per patient, gene ID mapping GDC v36 -> our v50, patient-level stratified hold-out (seed 42).
usage: python3 -I prepare_tcga.py <data_public/tcga_coad> <ref/gene_id_to_symbol.tsv> <out prep dir>   (stops with exit 2 on any PREREG stop rule)"""
import json, os, random, sys
import pandas as pd

tdir, sym_file, out = sys.argv[1:4]
if os.listdir(out): sys.exit(f"{out} is not empty: refusing to overwrite")
MAX_UNMAPPED_PCT = 5.0                                   # conf/thresholds.config analysis_max_unmapped_gene_pct (PREREG c)

def stop(msg): print("STOP (D1 FAIL):", msg); sys.exit(2)

# ---- (a) file selection
hits = json.load(open(os.path.join(tdir, "gdc_files_metadata.json")))["data"]["hits"]
rows = []
for h in hits:
    if len(h["cases"]) != 1 or len(h["cases"][0]["samples"]) != 1: stop(f"file {h['file_id']} has not exactly one case/sample")
    s = h["cases"][0]["samples"][0]
    grp = {"Primary Tumor": "tumour", "Solid Tissue Normal": "normal"}.get(s["sample_type"])
    if grp is None: stop(f"unexpected sample type {s['sample_type']}")
    rows.append(dict(file_id=h["file_id"], file_name=h["file_name"], patient_id=h["cases"][0]["submitter_id"], sample_id=s["submitter_id"], group=grp))
f = pd.DataFrame(rows).sort_values(["patient_id", "group", "sample_id", "file_id"])   # ASCII sort within patient and group
keep = f.groupby(["patient_id", "group"], sort=False).head(1)
drop = f.drop(keep.index).assign(reason="PREREG (a): patient has more than one file of this group; kept the first by (sample submitter_id, file_id)")
drop.to_csv(os.path.join(out, "dropped_files.tsv"), sep="\t", index=False)
if keep.sample_id.duplicated().any(): stop("two kept files share a sample ID")
keep = keep.sort_values("sample_id").reset_index(drop=True)
print(f"files {len(f)}; kept {len(keep)} ({keep.group.value_counts().to_dict()}); dropped {len(drop)}; patients {keep.patient_id.nunique()}")

# ---- (c) counts + mapping
ref = pd.read_csv(sym_file, sep="\t")
ref_u = ref.gene_id.str.split(".").str[0]
if ref_u.duplicated().any(): stop("our annotation has duplicated unversioned IDs")
u2v = dict(zip(ref_u, ref.gene_id))
cols, genes0 = {}, None
for _, r in keep.iterrows():
    d = pd.read_csv(os.path.join(tdir, "star_counts", r.file_name), sep="\t", comment="#", usecols=["gene_id", "gene_name", "unstranded"])
    d = d[~d.gene_id.str.startswith("N_") & ~d.gene_id.str.endswith("_PAR_Y")]
    if genes0 is None: genes0 = d[["gene_id", "gene_name"]].reset_index(drop=True)
    elif not d.gene_id.reset_index(drop=True).equals(genes0.gene_id): stop(f"gene rows differ in {r.file_name}")
    cols[r.sample_id] = d.unstranded.to_numpy()
cn = pd.DataFrame(cols, index=genes0.gene_id)
gu = genes0.gene_id.str.split(".").str[0]
mapped = gu.map(u2v)
unm = genes0[mapped.isna().to_numpy()]
pct = 100 * len(unm) / len(genes0)
cnt_pct = 100 * cn[mapped.isna().to_numpy()].to_numpy().sum() / cn.to_numpy().sum()
unm.to_csv(os.path.join(out, "unmapped_genes.tsv"), sep="\t", index=False)
summ = {"gdc_gene_rows_after_N_and_PAR_Y_removal": len(genes0), "unmapped_rows": len(unm), "unmapped_pct": round(pct, 4),
        "max_unmapped_pct": MAX_UNMAPPED_PCT, "unmapped_count_pct_of_total": round(cnt_pct, 4)}
print("mapping:", summ)
if pct > MAX_UNMAPPED_PCT: json.dump(summ, open(os.path.join(out, "mapping_summary.json"), "w"), indent=1); stop(f"unmapped {pct:.2f} % > {MAX_UNMAPPED_PCT} %")
cn = cn[mapped.notna().to_numpy()]; cn.index = mapped[mapped.notna()].to_numpy(); cn.index.name = "gene_id"
if cn.index.duplicated().any(): stop(f"{int(cn.index.duplicated().sum())} GDC rows map to the same v50 ID")
summ["genes_in_matrix"] = len(cn)
json.dump(summ, open(os.path.join(out, "mapping_summary.json"), "w"), indent=1)
cn.to_csv(os.path.join(out, "counts_tcga_coad_v50.tsv"), sep="\t")

# ---- (b) patient-level stratified split, seed 42
pg = keep.groupby("patient_id").group.agg(lambda g: ",".join(sorted(set(g))))
strat = pg.map({"normal,tumour": "paired", "tumour": "tumour_only", "normal": "normal_only"})
rng = random.Random(42); hold_p = set(); n_by = {}
for st in ["paired", "tumour_only", "normal_only"]:
    p = sorted(strat.index[strat == st]); rng.shuffle(p); n = round(len(p) / 3); hold_p |= set(p[:n]); n_by[st] = (len(p), n)
keep["stratum"] = keep.patient_id.map(strat); keep["split"] = ["holdout" if p in hold_p else "training" for p in keep.patient_id]
keep.to_csv(os.path.join(out, "split.tsv"), sep="\t", index=False)
open(os.path.join(out, "holdout_samples.txt"), "w").write("\n".join(keep.sample_id[keep.split == "holdout"]) + "\n")
cl = keep.assign(baseline=(keep.group == "normal").map({True: "TRUE", False: "FALSE"}), label=keep.sample_id)[["sample_id", "patient_id", "group", "baseline", "label"]]
cl.to_csv(os.path.join(out, "clinical_tcga_coad.tsv"), sep="\t", index=False)
print("strata (patients, hold-out patients):", n_by)
print(pd.crosstab(keep.split, keep.group))
