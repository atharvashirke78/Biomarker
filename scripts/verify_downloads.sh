#!/usr/bin/env bash
# Phase A verification: 40 files, gzip -t on each, R1/R2 read counts equal, table of run/group/sizes/reads.
set -uo pipefail
B=/workspace/Biomarker_pipeline
OUT=$B/logs/verify_downloads.tsv
printf "run\tgroup\tR1_bytes\tR2_bytes\tR1_reads\tR2_reads\tstatus\n" > $OUT
n=0; bad=0
while IFS=, read -r sample run study group; do
  [ "$sample" = sample ] && continue
  D=$B/data/$([ "$group" = tumour ] && echo tumor_sample || echo normal_sample)
  r1=$D/${run}_R1.fastq.gz; r2=$D/${run}_R2.fastq.gz
  st=OK
  for f in $r1 $r2; do
    if [ -s "$f" ]; then n=$((n+1)); gzip -t "$f" || st=GZIP_FAIL; else st=MISSING; fi
  done
  c1=NA; c2=NA
  if [ $st = OK ]; then
    c1=$(( $(zcat $r1 | wc -l) / 4 )); c2=$(( $(zcat $r2 | wc -l) / 4 ))
    [ "$c1" = "$c2" ] || st=READCOUNT_MISMATCH
  fi
  [ $st = OK ] || bad=$((bad+1))
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$run" "$group" "$(stat -c %s $r1 2>/dev/null || echo NA)" "$(stat -c %s $r2 2>/dev/null || echo NA)" "$c1" "$c2" "$st" >> $OUT
done < $B/samplesheets/diagnostic_20.csv
echo "files present: $n / 40 ; runs with problems: $bad"
column -t -s $'\t' $OUT
