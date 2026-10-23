#!/usr/bin/env python3
# Per-sample QC row. Written from the file formats of STAR Log.final.out, featureCounts .summary, Salmon meta_info.json, fastp json.
# usage: sample_qc.py id study group star_log fc_summary quant_dir fastp_json star_warn fc_warn salmon_warn star_fail salmon_fail
import json, os, re, sys
sid, study, group, star_log, fc_sum, quant, fjson = sys.argv[1:8]
star_warn, fc_warn, sal_warn, star_fail, sal_fail = map(float, sys.argv[8:13])

def star_unique(p):
    if os.path.basename(p).startswith("NO_FILE"): return None
    for l in open(p):
        if "Uniquely mapped reads %" in l:
            return float(l.split("|")[1].strip().rstrip("%"))
    return None
def fc_assigned(p):
    if os.path.basename(p).startswith("NO_FILE"): return None
    tot = ass = 0
    for i, l in enumerate(open(p)):
        if i == 0: continue
        k, v = l.rstrip("\n").split("\t")[:2]
        tot += int(v); ass += int(v) if k == "Assigned" else 0
    return 100.0 * ass / tot if tot else None
def salmon_mapped(q):
    return float(json.load(open(os.path.join(q, "aux_info", "meta_info.json")))["percent_mapped"])
def fastp(p):
    j = json.load(open(p))
    q30 = 100 * j["summary"]["after_filtering"]["q30_rate"]
    tot = j["summary"]["before_filtering"]["total_reads"]
    ok = j["filtering_result"]["passed_filter_reads"]
    return q30, 100.0 * ok / tot if tot else None

su, fa, sm = star_unique(star_log), fc_assigned(fc_sum), salmon_mapped(quant)
q30, passed = fastp(fjson)
status = "PASS"; reasons = []
def chk(name, v, lim, level):
    if v is not None and v < lim: reasons.append((level, "%s %.2f%% < %s limit %g%%" % (name, v, level, lim)))
chk("STAR unique", su, star_warn, "WARN"); chk("featureCounts assigned", fa, fc_warn, "WARN"); chk("Salmon mapped", sm, sal_warn, "WARN")
chk("STAR unique", su, star_fail, "FAIL"); chk("Salmon mapped", sm, sal_fail, "FAIL")
if reasons: status = "FAIL" if any(l == "FAIL" for l, _ in reasons) else "WARN"
reason = "; ".join(t for l, t in reasons if l == status) if reasons else ""
f = lambda x: "NA" if x is None else "%.2f" % x
# columns: id study group star_unique fc_assigned salmon_mapped fastp_q30 fastp_pass status reason (reason = which metric crossed which limit; empty for PASS)
print("\t".join([sid, study, group, f(su), f(fa), f(sm), f(q30), f(passed), status, reason]))
