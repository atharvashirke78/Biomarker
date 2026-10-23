#!/usr/bin/env bash
# QC RULE check: every results folder 01..08 must have qc/ with qc/fastp/<SRR>.fastp.{html,json} of ALL samples and that folder's multiqc_report.html.
# usage: check_qc_folders.sh <results_dir> <n_samples> [nonempty|exists]
R=${1:?results dir}; N=${2:?n samples}; MODE=${3:-nonempty}; bad=0
for d in 01_fastp 02_alignment 03_quant_salmon 04_counts 05_cohort_qc 06_dge 07_pathway 08_multiqc; do
  q=$R/$d/qc; h=$(ls $q/fastp/*.fastp.html 2>/dev/null | wc -l); j=$(ls $q/fastp/*.fastp.json 2>/dev/null | wc -l); m=$([ -e $q/multiqc_report.html ] && echo yes || echo no)
  ne=ok; if [ $MODE = nonempty ]; then [ -s $q/multiqc_report.html ] || ne=EMPTY_MULTIQC; [ $(find $q/fastp -name '*.fastp.*' -size 0 2>/dev/null | wc -l) -eq 0 ] || ne=EMPTY_FASTP; fi
  s=OK; { [ -d $q ] && [ $h -eq $N ] && [ $j -eq $N ] && [ $m = yes ] && [ $ne = ok ]; } || { s=PROBLEM; bad=$((bad+1)); }
  printf "%-16s qc/fastp/: html=%-3s json=%-3s multiqc_report.html=%-3s %s %s\n" $d $h $j $m $ne $s
done
[ $bad -eq 0 ] && echo "QC RULE: all 8 folders OK" || echo "QC RULE: $bad folder(s) with problems"; exit $bad
