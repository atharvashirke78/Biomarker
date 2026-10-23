#!/usr/bin/env bash
# Full download of all 20 runs, 3 at a time. Safe to rerun: finished files are skipped, partial files resume.
B=/workspace/Biomarker_pipeline
mkdir -p $B/logs
tail -n +2 $B/samplesheets/diagnostic_20.csv | awk -F, '{print $2, $4}' \
  | xargs -P 3 -L 1 bash $B/scripts/download_run.sh 2>&1 | tee -a $B/logs/download_full.log
