#!/usr/bin/env bash
# Per-step done counts (COMPLETED + CACHED, exit 0) from the newest trace file, tasks executed in this run with peak_rss,
# Salmon overlap check, tasks in flight, memory and disk.
# NOTE: a Nextflow trace file only gets a row when a task FINISHES, so "in flight" cannot come from the trace; it is
# derived as: submitted in the newest .nextflow.log but not yet in the trace.
cd /workspace/Biomarker_pipeline || exit 1
T=$(ls -t logs/nf_trace_*.tsv 2>/dev/null | head -1); [ -z "$T" ] && T=logs/nf_trace.tsv
echo "trace: $T  ($(date '+%F %T'))"
awk -F'\t' 'NR>1{split($4,a," \\("); p=a[1]; s=$5;
  if((s=="COMPLETED"||s=="CACHED")&&$6==0) d[p]++; else if(s=="FAILED"||s=="ABORTED") f[p]++; seen[p]=1}
  END{for(p in seen) printf "%-22s done=%d failed=%d\n", p, d[p], f[p]}' "$T" | sort
echo "executed in this run (not cached) - peak_rss:"
awk -F'\t' 'NR>1 && $5!="CACHED" && $4~/^(SALMON_QUANT|STAR_ALIGN)/{printf "  %-24s %-9s exit=%s submit=%s dur=%s peak_rss=%s\n",$4,$5,$6,$7,$8,$11}' "$T"
python3 -I - "$T" <<'PY'
import sys,csv,re,datetime as dt
def secs(s):
    t=0
    for n,u in re.findall(r'([\d.]+)\s*(ms|[dhms])',s): t+=float(n)*{'d':86400,'h':3600,'m':60,'s':1,'ms':.001}[u]
    return t
rows=[r for r in csv.DictReader(open(sys.argv[1]),delimiter='\t') if r['name'].startswith('SALMON_QUANT') and r['status']=='COMPLETED']
iv=[(dt.datetime.strptime(r['submit'][:19],'%Y-%m-%d %H:%M:%S').timestamp(),r['name']) for r in rows]
iv=[(a,a+secs(next(x['duration'] for x in rows if x['name']==n)),n) for a,n in iv]
ov=[(x[2],y[2]) for i,x in enumerate(iv) for y in iv[i+1:] if x[0]<y[1] and y[0]<x[1]]
print("Salmon overlap (executed this run):", ov if ov else "none" if rows else "no Salmon task finished yet")
PY
echo "in flight (from .nextflow.log: submitted, not in trace):"
grep -oE "Submitted process > [A-Z_]+ \([A-Za-z0-9_]+\)" .nextflow.log | sed 's/Submitted process > //' | sort -u > /tmp/.sub.$$ 2>/dev/null || sort -u </dev/null > /tmp/.sub.$$
awk -F'\t' 'NR>1{print $4}' "$T" | sort -u > /tmp/.done.$$
comm -23 /tmp/.sub.$$ /tmp/.done.$$ | sed 's/^/  /'; rm -f /tmp/.sub.$$ /tmp/.done.$$
free -g | sed -n 1,2p; df -h /workspace | tail -1
