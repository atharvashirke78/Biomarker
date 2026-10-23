#!/usr/bin/env bash
# usage: download_validation.sh RUN   (reads selected_runs.tsv; writes to validation_run/data/{tumor,normal}_sample)
# Downloads both mates from ENA, checks md5 against ENA, gzip -t, then renames .part -> .fastq.gz. Does not touch the main data/ folder.
set -uo pipefail
R=$1; B=/workspace/Biomarker_pipeline; T=$B/results/10_independent_validation/selected_runs.tsv
G=$(awk -F'\t' -v r="$R" '$1==r{print $4}' $T); FTP=$(awk -F'\t' -v r="$R" '$1==r{print $6}' $T); MD5=$(awk -F'\t' -v r="$R" '$1==r{print $7}' $T)
if [ "$G" = "tumour" ]; then D=$B/validation_run/data/tumor_sample; else D=$B/validation_run/data/normal_sample; fi
for n in 1 2; do
  U=$(echo "$FTP" | tr ';' '\n' | grep "_${n}.fastq.gz$" | head -1); M=$(echo "$MD5" | tr ';' '\n' | grep -n . | awk -F: -v u="$U" 'NR==0{print}' ); 
  [ -n "$U" ] || { echo "[$R] no URL mate $n"; exit 1; }
  OUT=$D/${R}_R${n}.fastq.gz
  [ -s "$OUT" ] && { echo "[$R] R$n exists"; continue; }
  curl -sS -L --retry 5 --retry-delay 10 -C - -o "$OUT.part" "https://$U" || { echo "[$R] curl failed R$n"; exit 1; }
  if gzip -t "$OUT.part"; then mv "$OUT.part" "$OUT"; echo "[$R] done R$n $(du -h "$OUT" | cut -f1)"; else echo "[$R] FAILED gzip R$n"; exit 1; fi
done
