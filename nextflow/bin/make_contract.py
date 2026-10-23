#!/usr/bin/env python3
# NEW (refactor Step 2, design section 2.3/2.4). Builds the contract sample_sheet.tsv and qc_table.tsv from the resolved QUANT sample sheet,
# the cohort sample_qc.tsv and (optional) cohort_qc_flags.tsv. FAIL samples are dropped from the sheet; rows are written in the order of the counts matrix header.
# usage: make_contract.py <resolved_samplesheet.csv (sample,run,study,group,batch,baseline)> <sample_qc.tsv> <cohort_qc_flags.tsv|NO_FILE> <matrix.tsv> <count_source> <out_sheet.tsv> <out_qc.tsv>
import csv, os, sys
sheet, sqc, flags, matrix, source, out_sheet, out_qc = sys.argv[1:8]
raw = {r["run"]: r for r in csv.DictReader(open(sheet, newline=""))}
qc = list(csv.DictReader(open(sqc, newline=""), delimiter="\t"))
fl = {} if os.path.basename(flags) == "NO_FILE" else {r["srr"]: r for r in csv.DictReader(open(flags, newline=""), delimiter="\t")}
mcols = open(matrix).readline().rstrip("\n").split("\t")[1:]
status = {r["srr"]: r["status"] for r in qc}
for s in mcols:
    if s not in raw: sys.exit("ERROR: matrix column %s not in the sample sheet" % s)
    if status.get(s) == "FAIL": sys.exit("ERROR: matrix column %s has QC status FAIL" % s)
with open(out_sheet, "w") as o:
    w = csv.writer(o, delimiter="\t", lineterminator="\n"); w.writerow(["sample_id", "label", "group", "baseline", "study", "batch", "count_source"])
    for s in mcols:
        r = raw[s]; w.writerow([s, r["sample"], r["group"], r["baseline"], r["study"], r["batch"], source])
with open(out_qc, "w") as o:
    w = csv.writer(o, delimiter="\t", lineterminator="\n")
    w.writerow(["sample_id", "label", "group", "study", "star_unique_pct", "fc_assigned_pct", "salmon_mapped_pct", "fastp_q30_pct", "fastp_pass_pct", "status", "reason", "flag_outlier_mad"])
    for r in sorted(qc, key=lambda r: r["srr"]):
        f = fl.get(r["srr"], {}).get("flag_outlier_3mad", "NA")
        w.writerow([r["srr"], r["label"], r["group"], r["study"], r["star_unique_pct"], r["fc_assigned_pct"], r["salmon_mapped_pct"], r["fastp_q30_pct"], r["fastp_pass_pct"], r["status"], r.get("reason", ""), f])
print("contract: %d samples in sheet, %d rows in qc_table (%d FAIL excluded from the sheet)" % (len(mcols), len(qc), sum(1 for r in qc if r["status"] == "FAIL")))
