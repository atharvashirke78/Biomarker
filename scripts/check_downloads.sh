#!/usr/bin/env bash
# For every run: read count in R1 and R2 (must be equal and match ENA's read_count), corruption check, and R1/R2 must not be
# copies of each other (first 10,000 reads; more than 50 % identical sequences -> PROBLEM; seen for SRR13401190).
# usage: check_downloads.sh [RUN ...]   (no argument: all 20 runs of the sample sheet; with arguments: only those runs)
#   env: CHECK_OUT=<file> to write the table somewhere else (default logs/integrity_check.tsv), SHEET=<csv sample,run,study,group>
B=/workspace/Biomarker_pipeline
OUT=${CHECK_OUT:-$B/logs/integrity_check.tsv}
SHEET=${SHEET:-$B/samplesheets/diagnostic_20.csv}

count() {
  f=$1
  [ -s "$f" ] || { echo MISSING; return; }
  t=$(mktemp)
  zcat "$f" | wc -l > "$t"; st=${PIPESTATUS[0]}
  if [ "$st" -ne 0 ]; then echo CORRUPT; else echo $(( $(cat "$t") / 4 )); fi
  rm -f "$t"
}

# prints "<identical> <total>" for the sequence lines of the first 10,000 reads of R1 and R2
ident() {
  paste <(zcat "$1" 2>/dev/null | head -n 40000 | awk 'NR%4==2') <(zcat "$2" 2>/dev/null | head -n 40000 | awk 'NR%4==2') \
    | awk -F'\t' '{ n++; if ($1 == $2) s++ } END { print s+0, n+0 }'
}

check_one() {
  R=$1; G=$2; B=/workspace/Biomarker_pipeline
  if [ "$G" = "tumour" ]; then D=$B/data/tumor_sample; else D=$B/data/normal_sample; fi
  EXP=$(awk -F'\t' -v r="$R" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $c["run_accession"]==r {print $c["read_count"]}' $B/samplesheets/ena_check.tsv)
  A=$(count "$D/${R}_R1.fastq.gz")
  C=$(count "$D/${R}_R2.fastq.gz")
  read SAME TOT < <(ident "$D/${R}_R1.fastq.gz" "$D/${R}_R2.fastq.gz")
  if [ "$A" = "$EXP" ] && [ "$C" = "$EXP" ] && [ "$TOT" -gt 0 ] && [ $(( 2 * SAME )) -le "$TOT" ]; then S=OK; else S=PROBLEM; fi
  printf "%s\t%s\t%s\t%s\t%s\t%s/%s\t%s\n" "$R" "$G" "$EXP" "$A" "$C" "$SAME" "$TOT" "$S"
}
export -f count check_one ident

printf "run\tgroup\tena_reads\tR1_reads\tR2_reads\tidentical_R1R2_first10k\tstatus\n" > "$OUT"
tail -n +2 "$SHEET" | awk -F, -v want="$*" 'BEGIN{n=split(want,w," "); for(i=1;i<=n;i++) k[w[i]]=1} n==0 || ($2 in k) {print $2, $4}' \
  | xargs -P 4 -L 1 bash -c 'check_one "$0" "$1"' >> "$OUT"
column -t "$OUT" 2>/dev/null || cat "$OUT"
echo "== runs with status OK: $(grep -c -P '\tOK$' "$OUT") of $(( $(wc -l < "$OUT") - 1 ))"
grep -q -P '\tPROBLEM$' "$OUT" && { echo "== PROBLEM runs: $(grep -P '\tPROBLEM$' "$OUT" | cut -f1 | tr '\n' ' ')"; exit 1; }
exit 0
