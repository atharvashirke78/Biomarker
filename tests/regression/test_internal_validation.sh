#!/usr/bin/env bash
# Tests of nextflow/bin/internal_validation.R (additive internal validation, upgrade H3). Read-only inputs, writes only into --work (a NEW folder).
# usage: tests/regression/test_internal_validation.sh --work <new folder> [--dev19 <run>/09_discover] [--demo <run>/09_discover] [--rscript <Rscript>]
# Cases: (1) the re-implemented recipe reproduces the pipeline model on the 19-sample data (genes, lambda.1se, cutoff, scores);
#        (2) the step runs (reduced repeats) and writes the three files with the expected columns and the limitation text;
#        (3) two runs give identical files; (4) NOT RUN rows: no model (demo), fewer than 10 samples, a class below 5; exit code always 0.
set -uo pipefail
HERE=$(cd "$(dirname "$0")/../.." && pwd); SCRIPT=$HERE/nextflow/bin/internal_validation.R
WORK= ; DEV=/workspace/run_dev19/Result/09_discover; DEMO=/workspace/biomarker_jobs/5d6ede1e1e43/output/09_discover; RS=/workspace/Biomarker_pipeline/env/prefix/r_classifier/bin/Rscript
while (( $# )); do case "$1" in --work) WORK=$2; shift 2;; --dev19) DEV=$2; shift 2;; --demo) DEMO=$2; shift 2;; --rscript) RS=$2; shift 2;; *) echo "unknown argument $1"; exit 2;; esac; done
[[ -n $WORK && ! -e $WORK ]] || { echo "usage: $0 --work <NEW folder>"; exit 2; }
mkdir -p "$WORK"; PASS=0; FAIL=0
check() { if eval "$2"; then echo "PASS  $1"; PASS=$((PASS+1)); else echo "FAIL  $1   [$2]"; FAIL=$((FAIL+1)); fi; }
iv() { local d=$1 o=$2; shift 2; "$RS" "$SCRIPT" --counts $d/00_input/counts_training.tsv --sheet $d/00_input/sample_sheet.discover.tsv --candidates $d/02_candidates/candidate_table.tsv \
  --genes $d/03_model/model_genes.tsv --spec $d/03_model/model_spec.txt --scores $d/03_model/development_scores.tsv --out "$o" "$@"; }
if [[ -d $DEV ]]; then
  iv $DEV $WORK/chk --check-only true > $WORK/chk.log 2>&1; rc=$?
  check "recipe reproduces the pipeline model (exit 0)" "[[ $rc == 0 ]]"
  check "genes, scores, cutoff, lambda.1se all TRUE" "grep -q 'genes TRUE, scores TRUE, cutoff TRUE.*lambda.1se TRUE' $WORK/chk.log"
  check "lambda.1se = 0.00957921662124062" "grep -q 'lambda.1se TRUE (0.00957921662124[0-9]* ' $WORK/chk.log"
  check "cutoff = 1.352325902" "grep -q 'cutoff TRUE (1.352325902 vs 1.352325902)' $WORK/chk.log"
  for t in a b; do iv $DEV $WORK/run_$t --repeats 3 --bootstrap 10 --subsamples 10 > $WORK/run_$t.log 2>&1; check "reduced run $t exit 0" "[[ $? == 0 ]]"; done
  for f in internal_validation_summary.tsv selection_frequency.tsv INTERNAL_VALIDATION.md; do check "$f written and identical in both runs" "[[ -s $WORK/run_a/$f ]] && cmp -s $WORK/run_a/$f $WORK/run_b/$f"; done
  check "summary header" "[[ \$(head -1 $WORK/run_a/internal_validation_summary.tsv) == \$'metric\tvalue\tci_low\tci_high\tlabel' ]]"
  check "frequency header" "[[ \$(head -1 $WORK/run_a/selection_frequency.tsv) == \$'gene_id\tsymbol\tdirection\tfrequency\tin_final_model' ]]"
  for m in apparent_AUC nestedCV_AUC nestedCV_sensitivity nestedCV_specificity bootstrap_corrected_AUC bootstrap_cutoff nestedCV_AUC_lambda_min; do check "metric $m present" "grep -q '^$m\s' $WORK/run_a/internal_validation_summary.tsv"; done
  check "all 9 final genes listed with in_final_model TRUE" "[[ \$(awk -F'\t' '\$5==\"TRUE\"' $WORK/run_a/selection_frequency.tsv | wc -l) == 9 ]]"
  check "limitation and 'NOT independent validation' in the text" "grep -q 'candidate screening stays OUTSIDE' $WORK/run_a/INTERNAL_VALIDATION.md && grep -q 'NOT independent validation' $WORK/run_a/INTERNAL_VALIDATION.md"
  # fewer than 10 samples / class below 5: a fitted-looking spec but a tiny sheet
  mkdir -p $WORK/tiny/00_input $WORK/tiny/02_candidates $WORK/tiny/03_model; cp $DEV/02_candidates/candidate_table.tsv $WORK/tiny/02_candidates/; cp $DEV/03_model/{model_genes.tsv,model_spec.txt,development_scores.tsv} $WORK/tiny/03_model/
  cut -f1-9 $DEV/00_input/counts_training.tsv > $WORK/tiny/00_input/counts_training.tsv; head -9 $DEV/00_input/sample_sheet.discover.tsv > $WORK/tiny/00_input/sample_sheet.discover.tsv
  iv $WORK/tiny $WORK/o_tiny > $WORK/tiny.log 2>&1; check "fewer than 10 samples: exit 0 and NOT RUN row" "[[ $? == 0 ]] && grep -q 'NOT RUN: only 8 samples' $WORK/o_tiny/internal_validation_summary.tsv"
  mkdir -p $WORK/cls; cp -r $WORK/tiny/02_candidates $WORK/tiny/03_model $WORK/cls/; mkdir -p $WORK/cls/00_input
  python3 - "$DEV" "$WORK/cls" <<'P'
import sys, csv
dev, out = sys.argv[1], sys.argv[2]
sh = list(csv.reader(open(dev + '/00_input/sample_sheet.discover.tsv'), delimiter='\t')); h = sh[0]; gi = h.index('group'); bi = h.index('baseline')
tum = [r for r in sh[1:] if r[bi].upper() != 'TRUE']; nor = [r for r in sh[1:] if r[bi].upper() == 'TRUE']; keep = nor + tum[:4]        # 4 tumours: below 5
open(out + '/00_input/sample_sheet.discover.tsv', 'w').write('\n'.join('\t'.join(r) for r in [h] + keep) + '\n')
cn = [l.rstrip('\n').split('\t') for l in open(dev + '/00_input/counts_training.tsv')]; ids = [r[0] for r in keep]; cols = [0] + [cn[0].index(i) for i in ids]
open(out + '/00_input/counts_training.tsv', 'w').write('\n'.join('\t'.join(r[c] for c in cols) for r in cn) + '\n')
P
  iv $WORK/cls $WORK/o_cls > $WORK/cls.log 2>&1; check "class below 5 samples: exit 0 and NOT RUN row" "[[ $? == 0 ]] && grep -q 'NOT RUN: a class has fewer than 5' $WORK/o_cls/internal_validation_summary.tsv"
else echo "SKIP  19-sample cases (no $DEV)"; fi
if [[ -d $DEMO ]]; then
  iv $DEMO $WORK/o_demo > $WORK/demo.log 2>&1; rc=$?
  check "demo (NOT FITTED): exit 0" "[[ $rc == 0 ]]"
  check "demo: 'NOT RUN: no model' in all three files" "grep -q 'NOT RUN: no model' $WORK/o_demo/internal_validation_summary.tsv && grep -q 'NOT RUN: no model' $WORK/o_demo/selection_frequency.tsv && grep -q 'NOT RUN: no model' $WORK/o_demo/INTERNAL_VALIDATION.md"
else echo "SKIP  demo cases (no $DEMO)"; fi
echo "internal validation tests: $PASS passed, $FAIL failed"; (( FAIL == 0 ))
