#!/usr/bin/env bash
# Read-only: for each of the 30 expected FASTQs report state (complete/part/missing), gzip -t and md5 vs ENA fastq_md5.
B=/workspace/Biomarker_pipeline; T=$B/results/10_independent_validation/selected_runs.tsv
one() {
  R=$1; n=$2; B=/workspace/Biomarker_pipeline; T=$B/results/10_independent_validation/selected_runs.tsv
  G=$(awk -F'\t' -v r="$R" '$1==r{print $4}' $T); MD=$(awk -F'\t' -v r="$R" '$1==r{print $7}' $T | cut -d';' -f$n); SZ=$(awk -F'\t' -v r="$R" '$1==r{print $8}' $T | cut -d';' -f$n)
  if [ "$G" = tumour ]; then D=$B/validation_run/data/tumor_sample; else D=$B/validation_run/data/normal_sample; fi
  F=$D/${R}_R${n}.fastq.gz
  if [ -s "$F" ]; then
    gzip -t "$F" 2>/dev/null && Z=gzip_ok || Z=gzip_FAIL
    M=$(md5sum "$F" | cut -d' ' -f1); [ "$M" = "$MD" ] && Mm=md5_ok || Mm=md5_DIFF
    printf "%s\tR%s\tcomplete\t%s\t%s\t%s\n" $R $n $(stat -c %s "$F") $Z $Mm
  elif [ -e "$F.part" ]; then printf "%s\tR%s\tdownloading(.part)\t%s of %s\t-\t-\n" $R $n $(stat -c %s "$F.part") $SZ
  else printf "%s\tR%s\tnot_started\t0 of %s\t-\t-\n" $R $n $SZ; fi
}
export -f one
for r in $(tail -n +2 $T | cut -f1); do echo "$r 1"; echo "$r 2"; done | xargs -P 4 -L 1 bash -c 'one $0 $1' | sort
