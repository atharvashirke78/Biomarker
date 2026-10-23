#!/usr/bin/env bash
# R-level test of the rule 3b refit plan (nextflow/bin/candidate_rules.R) on the 19-sample development design; no DESeq2 fit is run (the test hook --plan-only stops
# after robustness_refits.tsv and the settings line). The 19-sample clinical file has repeated patient_id values (paired tumour/normal), which is the case under test.
#   default                    -> the locked method: leave-one-SAMPLE-out, 19 refits, none skipped, loo_min_pass 18
#   --patient-level-cv true    -> leave-one-PATIENT-out: 11 refits (one per patient), pass >= ceil(0.9 x 11) = 10
# --plan-only is a TEST HOOK: it needs BIOMARKER_TEST_PLAN_ONLY=1 as well, is never passed by the workflow, and is checked below to appear nowhere in the pipeline code.
#
# usage: tests/integration/test_loo_refits.sh --work <NEW folder> [--install-dir /workspace/Biomarker_pipeline] [--inputs <discover out dir of a 19-sample run>]
# Every case prints PASS or FAIL; exit 0 only if all pass. Writes only into --work.
set -uo pipefail
HERE=$(cd "$(dirname "$0")/../.." && pwd)
WORK= ; INSTALL=${BIOMARKER_INSTALL_DIR:-/workspace/Biomarker_pipeline}
INPUTS=$INSTALL/results_rc1prep_19/regression_control_20261021/discover/out
while (( $# )); do case "$1" in --work) WORK=${2:-}; shift 2 ;; --install-dir) INSTALL=${2:-}; shift 2 ;; --inputs) INPUTS=${2:-}; shift 2 ;; *) echo "unknown argument $1"; exit 2 ;; esac; done
[[ -n $WORK && ! -e $WORK ]] || { echo "usage: $0 --work <new folder> [--install-dir DIR] [--inputs DIR]"; exit 2; }
RS=$INSTALL/env/prefix/r_dge/bin/Rscript; [[ -x $RS ]] || { echo "Rscript not found: $RS"; exit 2; }
for f in 00_input/counts_training.tsv 00_input/sample_sheet.discover.tsv 01_de/deseq2_results_all.tsv 01_de/deseq2_salmon_results_all.tsv; do [[ -e $INPUTS/$f ]] || { echo "missing input $INPUTS/$f"; exit 2; }; done
mkdir -p "$WORK"; WORK=$(cd "$WORK" && pwd)
PASS=0; FAIL=0
check() { if eval "$2"; then echo "PASS  $1"; PASS=$((PASS+1)); else echo "FAIL  $1   [$2]"; FAIL=$((FAIL+1)); fi; }
plan() {   # plan <folder> [extra args...]: runs the plan in <folder>, output in <folder>/out.txt, exit code in <folder>/rc
  local d=$WORK/$1; shift; mkdir -p "$d"
  (cd "$d" && BIOMARKER_TEST_PLAN_ONLY=${PLAN_ENV-1} "$RS" "$HERE/nextflow/bin/candidate_rules.R" --counts "$INPUTS/00_input/counts_training.tsv" --sheet "$INPUTS/00_input/sample_sheet.discover.tsv" \
    --de-primary "$INPUTS/01_de/deseq2_results_all.tsv" --de-secondary "$INPUTS/01_de/deseq2_salmon_results_all.tsv" --rules "$INSTALL/config/candidate_rules.yaml" \
    --thresholds "$INSTALL/conf/thresholds.config" --questions "$INSTALL/config/questions.yaml" "$@" > out.txt 2>&1; echo $? > rc)
}

# the premise: repeated patient_id values in the 19-sample design
check "premise: 19 samples, repeated patient_id (paired design)" "[[ \$(tail -n +2 $INPUTS/00_input/sample_sheet.discover.tsv | wc -l) == 19 && \$(tail -n +2 $INPUTS/00_input/sample_sheet.discover.tsv | cut -f8 | sort | uniq -d | wc -l) -gt 0 ]]"

# 1. default = the locked method
plan default --plan-only true
check "default: exit 0" "[[ \$(cat $WORK/default/rc) == 0 ]]"
check "default: 19 refits, all run (none skipped)" "[[ \$(tail -n +2 $WORK/default/robustness_refits.tsv | wc -l) == 19 && \$(tail -n +2 $WORK/default/robustness_refits.tsv | cut -f3 | sort -u) == run ]]"
check "default: refits are leave-one-SAMPLE-out (each refit = one SRR id, 18 samples)" "[[ \$(tail -n +2 $WORK/default/robustness_refits.tsv | cut -f1 | grep -c '^SRR') == 19 && \$(tail -n +2 $WORK/default/robustness_refits.tsv | cut -f2 | sort -u) == 18 ]]"
check "PASS default = locked method: 19 refits, loo_min_pass 18" "grep -q 'robustness loo, 19 refits, pass >= 18 (= ceil(0.9 x 19))' $WORK/default/out.txt && ! grep -q 'derived loo_min_pass' $WORK/default/out.txt"
check "default: repeated-patient message in the log, no patient-level wording" "grep -q 'Repeated patient_id found; folds and refits are sample-level (the locked method). Use --patient_level_cv for patient-level.' $WORK/default/out.txt && ! grep -qi 'leave-one-PATIENT' $WORK/default/out.txt"
plan explicit_false --patient-level-cv false --plan-only true
check "explicit false = default (identical robustness_refits.tsv)" "cmp -s $WORK/default/robustness_refits.tsv $WORK/explicit_false/robustness_refits.tsv"

# 2. the switch gives patient-level refits
plan patient --patient-level-cv true --plan-only true
check "switch true: exit 0" "[[ \$(cat $WORK/patient/rc) == 0 ]]"
check "switch true: 11 refits named by patient, pass >= 10" "[[ \$(tail -n +2 $WORK/patient/robustness_refits.tsv | wc -l) == 11 && \$(tail -n +2 $WORK/patient/robustness_refits.tsv | cut -f1 | grep -c '^[AB]_') == 11 ]] && grep -q 'robustness loo, 11 refits, pass >= 10' $WORK/patient/out.txt"
check "switch true: the locked-method guard (loo_min_pass 18) does not stop an opted-in patient-level plan" "! grep -q 'derived loo_min_pass' $WORK/patient/out.txt"

# 3. bad switch value, and --plan-only cannot be triggered by accident
plan badvalue --patient-level-cv maybe --plan-only true
check "switch value other than true/false is refused" "[[ \$(cat $WORK/badvalue/rc) != 0 ]]"
PLAN_ENV= plan noenv --plan-only true
check "--plan-only without BIOMARKER_TEST_PLAN_ONLY=1 is refused" "[[ \$(cat $WORK/noenv/rc) != 0 ]] && grep -q 'test hook' $WORK/noenv/out.txt"
check "--plan-only is in no workflow/module/config/report code" "! grep -rIl -e 'plan-only' -e 'plan_only' $HERE/modules $HERE/nextflow/main.nf $HERE/nextflow/nextflow.config $HERE/conf $HERE/config $HERE/scripts $HERE/nextflow/bin/package_model.py $HERE/nextflow/bin/fit_diagnostic_model.R $HERE/nextflow/bin/check_discover_input.py | grep -q ."
check "the switch defaults to false in nextflow.config" "grep -Eq '^\\s*patient_level_cv\\s*=\\s*false' $HERE/nextflow/nextflow.config"
check "module passes the switch to candidate_rules.R and fit_diagnostic_model.R" "[[ \$(grep -c -e '--patient-level-cv \${params.patient_level_cv}' $HERE/modules/discover.nf) == 2 ]]"
check "fold grouping in fit_diagnostic_model.R depends on the switch" "grep -q 'REPEATED_PATIENTS <- PLC == \"true\" && any(duplicated(pat))' $HERE/nextflow/bin/fit_diagnostic_model.R"

echo; echo "$PASS passed, $FAIL failed"; [[ $FAIL == 0 ]]
