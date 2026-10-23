#!/usr/bin/env bash
# Builds a small EXAMPLE input folder from DEVELOPMENT samples (never the validation samples): 3 normal + 3 tumour, the first 100,000 read pairs of each.
# It is for learning the workflow and testing its mechanics (a few minutes of compute per mode, about 100 MB); the cohort is far too small for biology.
#
# usage: examples/make_example_input.sh --out <NEW folder> [--pairs 100000] [--install-dir /workspace/Biomarker_pipeline]
# Then:  nextflow run .../nextflow/main.nf --input_dir <that folder> --mode apply --cohort development      (or the default discover mode)
# The 6 samples come from data/{normal_sample,tumor_sample} of the install (two public studies, both groups in both studies; all PASS the sample QC in the development run).
set -euo pipefail
OUT= ; PAIRS=100000 ; INSTALL=${BIOMARKER_INSTALL_DIR:-/workspace/Biomarker_pipeline}
while (( $# )); do case "$1" in --out) OUT=${2:-}; shift 2 ;; --pairs) PAIRS=${2:-}; shift 2 ;; --install-dir) INSTALL=${2:-}; shift 2 ;; *) echo "unknown argument $1"; exit 2 ;; esac; done
[[ -n $OUT && ! -e $OUT ]] || { echo "usage: $0 --out <new folder> [--pairs N] [--install-dir DIR]   (the folder must not exist)"; exit 2; }
# sample  group  study           (development samples; SRP301216 and SRP219837 each contribute normal and tumour)
SAMPLES=(
  "SRR13401193 normal SRP301216 normal_sample"
  "SRR13401194 normal SRP301216 normal_sample"
  "SRR10043507 normal SRP219837 normal_sample"
  "SRR13401203 tumour SRP301216 tumor_sample"
  "SRR10043499 tumour SRP219837 tumor_sample"
  "SRR10043500 tumour SRP219837 tumor_sample"
)
mkdir -p "$OUT"
echo 'sample_id,group,batch_or_study' > "$OUT/metadata.csv"
for row in "${SAMPLES[@]}"; do
  read -r s g st d <<< "$row"
  for m in R1 R2; do
    src="$INSTALL/data/$d/${s}_$m.fastq.gz"; [[ -s $src ]] || { echo "ERROR: $src not found"; exit 1; }
    set +o pipefail; zcat "$src" | head -n $((PAIRS * 4)) | gzip -1 > "$OUT/${s}_$m.fastq.gz"; set -o pipefail
  done
  echo "$s,$g,$st" >> "$OUT/metadata.csv"
  echo "  $s ($g, $st): $PAIRS read pairs"
done
echo "example input folder: $OUT ($(du -sh "$OUT" | cut -f1), 6 samples)"
