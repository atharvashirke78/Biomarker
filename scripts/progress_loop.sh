#!/usr/bin/env bash
# Append a status line to docs/run_log.md every 25 min while the full run is alive (survives tmux loss, not a pod stop).
cd /workspace/Biomarker_pipeline || exit 1
while pgrep -f "nextflow/main.nf" >/dev/null || pgrep -f "env/mamba/envs/nf" >/dev/null; do
  sleep 1500
  { echo "- $(date '+%F %H:%M') progress: $(bash scripts/status.sh | grep -E 'done=' | sed 's/  */ /g' | tr '\n' ';') | $(free -g | awk 'NR==2{print "mem used "$3" GB free "$4" GB"}') | $(df -h /workspace | awk 'NR==2{print "disk avail "$4" ("$5" used)"}')"; } >> docs/run_log.md
done
echo "- $(date '+%F %H:%M') run process no longer alive; see logs/nf_full_*.log" >> docs/run_log.md
