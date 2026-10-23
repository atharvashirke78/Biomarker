#!/usr/bin/env python3
# Builds the validation sample sheet + an exclusion table from the ENA filereports (metadata only).
# Selection: paper Table S2 runs (SRP344867: 5 cancer + 4 normal; SRP245232: 3 normal + 3 primary). Each ENA run gets a reason.
import csv, sys
B = "/workspace/Biomarker_pipeline/results/10_independent_validation"
S2 = {  # runs listed in paper Table S2 (from step 1 comparison)
 "SRP344867": {"SRR16832142","SRR16832141","SRR16832130","SRR16832119","SRR16832118","SRR16832116","SRR16832115","SRR16832114","SRR16832113"},
 "SRP245232": {"SRR10967902","SRR10967903","SRR10967905","SRR10967906","SRR10967908","SRR10967909"}}
sheet, excl = [], []
for st in S2:
    rows = list(csv.DictReader(open(f"{B}/{st}/metadata/ena_filereport.tsv"), delimiter="\t"))
    for r in rows:
        run, alias = r["run_accession"], r["sample_alias"]; t = (r["experiment_title"] + " " + r["tissue_type"]).lower()
        ok_type = r["library_strategy"] == "RNA-Seq" and r["library_source"] == "TRANSCRIPTOMIC" and r["library_layout"] == "PAIRED"
        if not ok_type: why = f"not paired-end mRNA/RNA-seq ({r['library_strategy']}, {r['library_layout']})"
        elif "adenoma" in t: why = "adenoma (not cancer or normal; not in paper S2)"
        elif "metast" in t: why = "liver metastasis (not in paper S2)"
        elif run not in S2[st]: why = "normal not in paper Table S2 (labelled sensitivity candidate)" if "normal" in t else "not in paper Table S2"
        else: why = None
        if why: excl.append((st, run, alias, r["experiment_title"], why)); continue
        if "normal" in t: g = "normal"
        elif "cancer" in t: g = "tumour"
        else: sys.exit(f"unclear label {run} {alias} {t}")
        sheet.append((run, run, st, g, r["read_count"], r["fastq_ftp"], r["fastq_md5"], r["fastq_bytes"], alias, r["experiment_title"]))
with open("/workspace/Biomarker_pipeline/validation_run/samplesheets/validation_15.csv", "w") as f:
    f.write("sample,run,study,group\n")
    for s in sheet: f.write(",".join(s[:4]) + "\n")
with open(f"{B}/selected_runs.tsv", "w") as f:
    f.write("run\tsample_alias\tstudy\tgroup\tena_read_count\tfastq_ftp\tfastq_md5\tfastq_bytes\texperiment_title\n")
    for s in sheet: f.write("\t".join([s[0], s[8], s[2], s[3], s[4], s[5], s[6], s[7], s[9]]) + "\n")
with open(f"{B}/excluded_runs.tsv", "w") as f:
    f.write("study\trun\tsample_alias\texperiment_title\treason\n")
    for e in excl: f.write("\t".join(e) + "\n")
from collections import Counter
print(Counter((s[2], s[3]) for s in sheet)); print(len(excl), "excluded")
