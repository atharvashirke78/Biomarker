#!/usr/bin/env python3
"""Download public counts-level data (no analysis). usage: python3 -I fetch_public.py <data_public dir>
TCGA-COAD (GDC open access): STAR - Counts for primary tumor + solid tissue normal, file metadata, case clinical/survival (/cases, expanded).
GSE78220 (GEO): NCBI-generated raw counts and the series matrix.
Every downloaded file -> one row in SOURCES.tsv (url, UTC date, size, sha256; GDC md5 checked). A failed download is logged and NOT retried."""
import datetime, hashlib, json, os, sys, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

out = os.path.abspath(sys.argv[1]); gdc = "https://api.gdc.cancer.gov"
tdir = os.path.join(out, "tcga_coad"); cdir = os.path.join(tdir, "star_counts"); gdir = os.path.join(out, "gse78220")
for d in (cdir, gdir): os.makedirs(d, exist_ok=True)
src = os.path.join(out, "SOURCES.tsv"); fail = os.path.join(out, "FAILED_DOWNLOADS.tsv")
cols = ["path", "url", "download_utc", "size_bytes", "sha256", "md5", "md5_expected", "md5_check", "note"]
if os.path.exists(src): sys.exit(f"{src} exists: refusing to overwrite")

def now(): return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def fetch(url, path, md5_expected="", note=""):
    if os.path.exists(path): return {"path": path, "url": url, "error": "target exists (not overwritten)"}
    t = now(); h256, h5, n = hashlib.sha256(), hashlib.md5(), 0
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "curl/8"}), timeout=300) as r, open(path + ".part", "wb") as f:
            for b in iter(lambda: r.read(1 << 20), b""): f.write(b); h256.update(b); h5.update(b); n += len(b)
    except Exception as e:
        if os.path.exists(path + ".part"): os.remove(path + ".part")
        return {"path": path, "url": url, "error": f"{type(e).__name__}: {e}"}
    os.rename(path + ".part", path)
    chk = "" if not md5_expected else ("OK" if h5.hexdigest() == md5_expected else "MISMATCH")
    return {"path": os.path.relpath(path, out), "url": url, "download_utc": t, "size_bytes": n, "sha256": h256.hexdigest(),
            "md5": h5.hexdigest(), "md5_expected": md5_expected, "md5_check": chk, "note": note}

def q(endpoint, params): return f"{gdc}/{endpoint}?" + urllib.parse.urlencode(params)

res = []
status = json.load(urllib.request.urlopen(f"{gdc}/status", timeout=60))
rel = f"GDC {status.get('data_release', '?')}"
fl = {"op": "and", "content": [
    {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-COAD"}},
    {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
    {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
    {"op": "=", "content": {"field": "access", "value": "open"}},
    {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor", "Solid Tissue Normal"]}}]}
ff = "file_id,file_name,md5sum,file_size,cases.case_id,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.tissue_type"
res.append(fetch(q("files", {"filters": json.dumps(fl), "fields": ff, "size": 2000, "format": "JSON"}),
                 os.path.join(tdir, "gdc_files_metadata.json"), note=f"{rel}; file list + md5 for STAR - Counts, primary tumor + solid tissue normal"))
meta = json.load(open(os.path.join(tdir, "gdc_files_metadata.json")))["data"]["hits"]
cl = {"op": "=", "content": {"field": "project.project_id", "value": "TCGA-COAD"}}
res.append(fetch(q("cases", {"filters": json.dumps(cl), "expand": "demographic,diagnoses,diagnoses.treatments,exposures,follow_ups",
                             "size": 2000, "format": "JSON"}),
                 os.path.join(tdir, "gdc_cases_clinical.json"), note=f"{rel}; all TCGA-COAD cases, demographic/diagnoses/follow_ups expanded (vital_status, days_to_death, days_to_last_follow_up)"))
for u, fn, nt in [
    ("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE78220&format=file&file=GSE78220_raw_counts_GRCh38.p13_NCBI.tsv.gz",
     "GSE78220_raw_counts_GRCh38.p13_NCBI.tsv.gz", "GEO NCBI-generated raw counts"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/matrix/GSE78220_series_matrix.txt.gz",
     "GSE78220_series_matrix.txt.gz", "GEO series matrix")]:
    res.append(fetch(u, os.path.join(gdir, fn), note=nt))

def one(h):
    s = h["cases"][0]["samples"][0]
    return fetch(f"{gdc}/data/{h['file_id']}", os.path.join(cdir, h["file_name"]), h["md5sum"],
                 f"{rel}; case {h['cases'][0]['submitter_id']}; sample {s['submitter_id']}; {s['sample_type']}; file_id {h['file_id']}")
with ThreadPoolExecutor(6) as ex: res += list(ex.map(one, meta))

ok = [r for r in res if "error" not in r]; bad = [r for r in res if "error" in r]
with open(src, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in ok: f.write("\t".join(str(r[c]) for c in cols) + "\n")
with open(fail, "w") as f:
    f.write("path\turl\terror\n")
    for r in bad: f.write(f"{r['path']}\t{r['url']}\t{r['error']}\n")
print(f"{rel}; files in GDC list: {len(meta)}; downloaded OK: {len(ok)}; failed: {len(bad)}; md5 mismatch: {sum(r['md5_check'] == 'MISMATCH' for r in ok)}")
