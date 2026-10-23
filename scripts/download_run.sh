#!/usr/bin/env bash
# usage: download_run.sh SRR_NUMBER GROUP   (GROUP = tumour or normal)
set -uo pipefail
R=$1; G=$2
B=/workspace/Biomarker_pipeline
if [ "$G" = "tumour" ]; then D=$B/data/tumor_sample; else D=$B/data/normal_sample; fi
mkdir -p "$D"
FTP=$(awk -F'\t' -v r="$R" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["run_accession"]==r {print $c["fastq_ftp"]}' $B/samplesheets/ena_check.tsv)
for n in 1 2; do
  U=$(echo "$FTP" | tr ';' '\n' | grep "_${n}.fastq.gz$" | head -1)
  if [ -z "$U" ]; then echo "[$R] no URL for mate $n"; exit 1; fi
  OUT=$D/${R}_R${n}.fastq.gz
  if [ -s "$OUT" ]; then echo "[$R] R$n exists, skipping"; continue; fi
  curl -sS -L --retry 5 --retry-delay 10 -C - -o "$OUT.part" "https://$U"
  if gzip -t "$OUT.part"; then mv "$OUT.part" "$OUT"; echo "[$R] done R$n $(du -h "$OUT" | cut -f1)"; else echo "[$R] FAILED R$n"; fi
done
