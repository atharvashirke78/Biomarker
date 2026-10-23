#!/usr/bin/env bash
# For every selected run: md5 of both mates vs ENA fastq_md5, read count R1/R2 vs ENA read_count, and R1/R2 identical-sequence test (first 10,000 reads;
# >50 % identical -> PROBLEM, as for SRR13401190). Output: results/10_independent_validation/integrity_check.tsv
B=/workspace/Biomarker_pipeline; T=$B/results/10_independent_validation/selected_runs.tsv; OUT=$B/results/10_independent_validation/integrity_check.tsv
one() {
  R=$1; B=/workspace/Biomarker_pipeline; T=$B/results/10_independent_validation/selected_runs.tsv
  G=$(awk -F'\t' -v r="$R" '$1==r{print $4}' $T); EXP=$(awk -F'\t' -v r="$R" '$1==r{print $5}' $T); MD=$(awk -F'\t' -v r="$R" '$1==r{print $7}' $T)
  if [ "$G" = "tumour" ]; then D=$B/validation_run/data/tumor_sample; else D=$B/validation_run/data/normal_sample; fi
  M1=$(md5sum $D/${R}_R1.fastq.gz | cut -d' ' -f1); M2=$(md5sum $D/${R}_R2.fastq.gz | cut -d' ' -f1)
  E1=$(echo $MD | cut -d';' -f1); E2=$(echo $MD | cut -d';' -f2)
  N1=$(( $(zcat $D/${R}_R1.fastq.gz | wc -l) / 4 )); N2=$(( $(zcat $D/${R}_R2.fastq.gz | wc -l) / 4 ))
  read SAME TOT < <(paste <(zcat $D/${R}_R1.fastq.gz | head -n 40000 | awk 'NR%4==2') <(zcat $D/${R}_R2.fastq.gz | head -n 40000 | awk 'NR%4==2') | awk -F'\t' '{n++; if($1==$2)s++} END{print s+0, n+0}')
  S=OK; [ "$M1" = "$E1" ] && [ "$M2" = "$E2" ] && [ "$N1" = "$EXP" ] && [ "$N2" = "$EXP" ] && [ "$TOT" -gt 0 ] && [ $((2*SAME)) -le "$TOT" ] || S=PROBLEM
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\t%s\n" $R $G $EXP $N1 $N2 "$([ $M1 = $E1 ] && [ $M2 = $E2 ] && echo md5_ok || echo md5_DIFF)" "" $SAME $TOT $S
}
export -f one
printf "run\tgroup\tena_reads\tR1_reads\tR2_reads\tmd5\t-\tidentical_R1R2_first10k\tstatus\n" > $OUT
tail -n +2 $T | cut -f1 | xargs -P 4 -n 1 bash -c 'one $0' >> $OUT
cat $OUT
