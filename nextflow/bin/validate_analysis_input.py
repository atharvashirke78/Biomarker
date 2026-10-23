#!/usr/bin/env python3
# NEW (refactor Step 2, docs/two_pipeline_design.md section 3). ANALYSIS input checks; stdlib only.
# Every failed check prints one line "CHECK <n> <name>: found ...; required ..." and the script exits 2 (nothing partial is published).
# usage: validate_analysis_input.py --counts counts.tsv --sheet sample_sheet.tsv --symbols gene_id_to_symbol.tsv --route featurecounts|salmon_tximport
#        --min_samples N --max_unmapped_pct P [--qc_table qc_table.tsv] [--secondary true|false] [--geneset_md5 file] --out_sheet sample_sheet.validated.tsv --out_json check.json
import argparse, csv, hashlib, json, math, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True); ap.add_argument("--sheet", required=True); ap.add_argument("--symbols", required=True)
ap.add_argument("--route", required=True); ap.add_argument("--min_samples", type=int, required=True); ap.add_argument("--max_unmapped_pct", type=float, required=True)
ap.add_argument("--qc_table", default=""); ap.add_argument("--secondary", default="false"); ap.add_argument("--geneset_md5", default="")
ap.add_argument("--out_sheet", required=True); ap.add_argument("--out_json", required=True)
a = ap.parse_args()
if os.path.realpath(a.out_sheet) == os.path.realpath(a.sheet): sys.exit('ERROR: --out_sheet must not be the input sheet')

ROUTES = ("featurecounts", "salmon_tximport")
fails, results = [], {}
def fail(n, name, found, required):
    fails.append("CHECK %d %s: found %s; required %s" % (n, name, found, required)); results[n] = {"name": name, "status": "FAIL", "found": found}
def ok(n, name, info=""):
    results.setdefault(n, {"name": name, "status": "OK", "info": info})
def lst(x, k=10):
    x = sorted(x); return ", ".join(x[:k]) + (" ... (+%d more)" % (len(x) - k) if len(x) > k else "")

# ---- read the sample sheet
with open(a.sheet, newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
need = ["sample_id", "group", "baseline", "batch", "count_source"]
cols = list(rows[0].keys()) if rows else []
if not rows or [c for c in need if c not in cols]:
    sys.exit("CHECK 4 sheet_columns: found columns [%s]; required %s" % (", ".join(cols), ", ".join(need)))
ids = [r["sample_id"] for r in rows]

# ---- 1 baseline group
bl = {}
for r in rows: bl.setdefault(r["group"], set()).add(r["baseline"].strip().upper())
bad = [g for g, v in bl.items() if not v <= {"TRUE", "FALSE"} or len(v) > 1]
base_groups = [g for g, v in bl.items() if v == {"TRUE"}]
if bad: fail(1, "baseline_flag", "group(s) with missing/inconsistent baseline values: " + lst(bad), "baseline TRUE or FALSE, identical for all rows of a group")
elif len(base_groups) == 0: fail(1, "baseline_group_missing", "no group with baseline = TRUE (groups: %s)" % lst(bl), "exactly one baseline group")
elif len(base_groups) > 1: fail(1, "baseline_group_ambiguous", "baseline = TRUE for %s" % lst(base_groups), "exactly one baseline group")
else: ok(1, "baseline_group", base_groups[0])

# ---- 7 group names come from the sheet; v1 supports exactly two groups
if len(bl) != 2: fail(7, "number_of_groups", "%d group(s): %s" % (len(bl), lst(bl)), "exactly 2 groups (baseline + one comparison group; more than two is not supported in v1)")
else: ok(7, "groups_from_sheet", lst(bl))

# ---- 4 sample sheet vs matrix
with open(a.counts, newline="") as fh:
    rd = csv.reader(fh, delimiter="\t"); header = next(rd)
    mcols = header[1:]
    genes, nonint, neg, nas, nrows = [], 0, 0, 0, 0
    first_bad = None
    for line in rd:
        genes.append(line[0]); nrows += 1
        for v in line[1:]:
            try: x = float(v)
            except ValueError:
                nas += 1; continue
            if math.isnan(x) or math.isinf(x): nas += 1; continue
            if x < 0: neg += 1
            if a.route == "featurecounts" and x != int(x):
                nonint += 1
                if first_bad is None: first_bad = (line[0], v)
dup_ids = sorted({i for i in ids if ids.count(i) > 1}); dup_cols = sorted({c for c in mcols if mcols.count(c) > 1})
only_sheet = sorted(set(ids) - set(mcols)); only_matrix = sorted(set(mcols) - set(ids))
if dup_ids or dup_cols or only_sheet or only_matrix:
    m = []
    if dup_ids: m.append("duplicated sample_id in sheet: " + lst(dup_ids))
    if dup_cols: m.append("duplicated columns in matrix: " + lst(dup_cols))
    if only_sheet: m.append("in sheet but not in matrix: " + lst(only_sheet))
    if only_matrix: m.append("in matrix but not in sheet: " + lst(only_matrix))
    fail(4, "sheet_matrix_mismatch", "; ".join(m), "the same sample IDs in the sample sheet and the matrix header")
elif ids != mcols:
    fail(4, "sheet_matrix_order", "sheet order and matrix column order differ", "same order (sheet rows listed in matrix column order)")
else: ok(4, "sheet_matrix_match", "%d samples" % len(ids))
if a.qc_table:
    with open(a.qc_table, newline="") as fh: q = {r.get("srr") or r.get("sample_id"): r for r in csv.DictReader(fh, delimiter="\t")}
    failed_in = sorted(i for i in ids if q.get(i, {}).get("status") == "FAIL")
    if failed_in: fail(4, "qc_fail_in_matrix", "sample(s) with QC status FAIL still in the analysis: " + lst(failed_in), "FAIL samples excluded from matrix and sheet")

# ---- 2 minimum samples per group
cnt = {g: sum(1 for r in rows if r["group"] == g) for g in bl}
low = {g: n for g, n in cnt.items() if n < a.min_samples}
if low: fail(2, "too_few_samples", ", ".join("group '%s' has %d" % (g, n) for g, n in sorted(low.items())), ">= %d samples per group (analysis_min_samples_per_group)" % a.min_samples)
else: ok(2, "min_samples_per_group", json.dumps(cnt))

# ---- 3 counts content (integer check only for count_source = featurecounts)
dup_genes = sorted({g for g in genes if genes.count(g) > 1}) if len(genes) != len(set(genes)) else []
msgs = []
if nonint: msgs.append("%d non-integer value(s) in a featurecounts matrix, first: gene %s value %s" % (nonint, first_bad[0], first_bad[1]))
if neg: msgs.append("%d negative value(s)" % neg)
if nas: msgs.append("%d NA/non-numeric value(s)" % nas)
if dup_genes: msgs.append("duplicate gene IDs: " + lst(dup_genes))
if msgs: fail(3, "counts_content", "; ".join(msgs), "%sno negatives, NAs or duplicate gene IDs" % ("integers (count_source = featurecounts), " if a.route == "featurecounts" else "(non-integer tximport counts are allowed) "))
else: ok(3, "counts_content", "%d genes x %d samples, route %s%s" % (nrows, len(mcols), a.route, ", integer check applied" if a.route == "featurecounts" else ", integer check not applied"))

# ---- 5 count_source recorded, known, uniform
src = sorted({r["count_source"].strip() for r in rows})
if "" in src: fail(5, "count_source_missing", "empty count_source for %d sample(s)" % sum(1 for r in rows if not r["count_source"].strip()), "one of %s for every sample" % "|".join(ROUTES))
elif [s for s in src if s not in ROUTES]: fail(5, "count_source_unknown", "value(s) " + lst([s for s in src if s not in ROUTES]), "one of %s" % "|".join(ROUTES))
elif len(src) > 1: fail(5, "count_source_mixed", "mix of " + lst(src), "one count source per analysis")
elif a.route not in ROUTES: fail(5, "route_unknown", a.route, "one of %s" % "|".join(ROUTES))
elif src[0] != a.route and a.secondary.lower() != "true": fail(5, "count_source_vs_route", "sheet says %s, analysis route is %s" % (src[0], a.route), "equal (the sheet records how the counts were produced)")
else: ok(5, "count_source", "sheet=%s route=%s%s" % (src[0], a.route, " (secondary cross-check route)" if a.secondary.lower() == "true" else ""))

# ---- 6 design estimable (H2b gate): no empty batch x group cell (= no batch nested in a group), and the DESeq2 model matrix (~ batch + group, or ~ group with one batch)
# has full column rank with residual degrees of freedom left. Only turns the error DESeq2 would raise later into an early, explained stop.
def model_rank(rws):
    bs = sorted({r["batch"] for r in rws}); gs = sorted({r["group"] for r in rws}, key=lambda g: (g not in base_groups, g))
    m = [[1.0] + [1.0 if r["batch"] == b else 0.0 for b in bs[1:]] + [1.0 if r["group"] == g else 0.0 for g in gs[1:]] for r in rws]
    ncol = len(m[0]); rank = 0; m = [row[:] for row in m]
    for c in range(ncol):
        piv = next((i for i in range(rank, len(m)) if abs(m[i][c]) > 1e-9), None)
        if piv is None: continue
        m[rank], m[piv] = m[piv], m[rank]
        for i in range(len(m)):
            if i != rank and abs(m[i][c]) > 1e-12:
                f = m[i][c] / m[rank][c]; m[i] = [x - f * y for x, y in zip(m[i], m[rank])]
        rank += 1
    return rank, ncol
cell = {}
for r in rows: cell[(r["batch"], r["group"])] = cell.get((r["batch"], r["group"]), 0) + 1
batches = sorted({r["batch"] for r in rows}); empty = [(b, g) for b in batches for g in bl if (b, g) not in cell]
gtab = sorted(bl)
table_txt = "batch-by-group table [%s]: " % " | ".join(gtab) + "; ".join("%s: %s" % (b, " | ".join(str(cell.get((b, g), 0)) for g in gtab)) for b in batches)
if empty: fail(6, "design_not_estimable", "empty batch x group cell(s): " + "; ".join("batch %s / group %s" % e for e in empty) + ". " + table_txt + ". The design ~ batch + group cannot be fitted when a batch holds only one group (group is then nested in batch)",
               "every batch has samples in both groups so ~ batch + group is estimable. To change: add samples of the missing group to that batch, merge batches (same batch label) if they are the same experiment, or remove the batch")
else:
    rk, nc = model_rank(rows)
    if rk < nc or len(rows) - rk < 1:
        fail(6, "design_not_estimable", "%s; the model matrix has rank %d of %d columns with %d samples (%d residual degrees of freedom)" % (table_txt, rk, nc, len(rows), len(rows) - rk),
             "a full-rank model matrix with at least 1 residual degree of freedom. To change: add samples, or merge batches so there are fewer batch terms than samples allow")
    else: ok(6, "design_estimable", "%d batch(es), cells: %s" % (len(batches), json.dumps({"%s|%s" % k: v for k, v in sorted(cell.items())})))

# ---- 9 gene IDs found in the gene-id-to-symbol table
with open(a.symbols, newline="") as fh:
    rd = csv.reader(fh, delimiter="\t"); next(rd); sym = {l[0] for l in rd}
miss = [g for g in genes if g not in sym]; pct = 100.0 * len(miss) / max(1, len(genes))
if pct > a.max_unmapped_pct: fail(9, "gene_ids_unmapped", "%.2f%% of %d gene IDs not in the symbol table (e.g. %s)" % (pct, len(genes), lst(miss, 3)), "<= %g%% unmapped (analysis_max_unmapped_gene_pct)" % a.max_unmapped_pct)
else: ok(9, "gene_ids_mapped", "%.2f%% unmapped (%d of %d)" % (pct, len(miss), len(genes)))

# ---- 8 gene-set files present with the expected md5 (file format: "<md5>  <path>" per line; skipped when no file is given)
if a.geneset_md5:
    bad = []
    for l in open(a.geneset_md5):
        if not l.strip(): continue
        md5, path = l.split(None, 1); path = path.strip()
        if not os.path.exists(path): bad.append("missing " + path); continue
        h = hashlib.md5(open(path, "rb").read()).hexdigest()
        if h != md5: bad.append("md5 differs " + path)
    if bad: fail(8, "geneset_files", "; ".join(bad), "all gene-set files present with the md5 recorded in the config (no silent fallback, no download)")
    else: ok(8, "geneset_files", "all present, md5 match")
else: ok(8, "geneset_files", "not requested for this call")

json.dump({"route": a.route, "secondary": a.secondary, "n_samples": len(ids), "baseline_group": base_groups[0] if len(base_groups) == 1 else None,
           "groups": cnt, "min_samples_per_group": a.min_samples, "checks": [results[k] for k in sorted(results)], "passed": not fails}, open(a.out_json, "w"), indent=2)
if fails:
    sys.stderr.write("ANALYSIS INPUT CHECK FAILED (%d):\n%s\n" % (len(fails), "\n".join(fails))); sys.exit(2)
with open(a.sheet) as src_f, open(a.out_sheet, "w") as o: o.write(src_f.read())
print("ANALYSIS input checks passed: %d samples, route %s, groups %s" % (len(ids), a.route, json.dumps(cnt)))
