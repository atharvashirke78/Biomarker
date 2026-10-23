#!/usr/bin/env bash
# Resume the 20-sample diagnostic full run. Report/timeline/trace/console log go to TIMESTAMPED files in logs/ (never overwritten).
cd /workspace/Biomarker_pipeline || exit 1
SID=$(cat docs/full_run_session_id.txt)
TS=$(date +%Y%m%d_%H%M%S)
echo "$TS" > logs/latest_run_ts.txt
env/bin/nextflow run nextflow/main.nf -profile conda -resume "$SID" \
  -with-report   logs/nf_report_$TS.html \
  -with-timeline logs/nf_timeline_$TS.html \
  -with-trace    logs/nf_trace_$TS.tsv \
  "$@" 2>&1 | tee logs/nf_full_$TS.log
