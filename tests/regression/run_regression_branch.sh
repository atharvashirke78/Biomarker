#!/usr/bin/env bash
# Two-level regression (owner decision 2026-10-16), part 1: run the reference workloads into a NEW folder.
#   A  discover/    -entry DISCOVER, 19 development samples, both routes, leave-one-out   (command: docs/run_log.md, 2026-10-15 Phase 1)
#   B  diagnostic/  scripts/run_diagnostic.sh  (ANALYSIS, locked model on development data, rule 4 self-consistency, 48 checks)
#   C  stub_default/  default entry (FASTQ sheet -> QUANT -> ANALYSIS) with -stub-run, own ref/work dirs (structure of the QUANT path)
#   extra workloads are added by later steps (stub runs of the folder entry), each in its own subfolder.
# Then compare with tests/regression/check_pod_baseline.py (byte level vs the pod baseline) and tests/regression/check_decisions.py
# (decision level vs the locked results, with an independent recomputation of the locked scores).
# usage: tests/regression/run_regression_branch.sh <new dir>
# BRANCH COPY of tests/regression/run_regression.sh (the original is unchanged). CODE = the checkout under test (default: the checkout this script is in);
# DATA = the install with results, environments and reference (default /workspace/Biomarker_pipeline). The ONLY other difference: workload A is
# started with --questions_config <DATA>/config/questions.yaml (byte-identical to CODE/config/questions.yaml, sha256 checked below), because the
# baseline provenance file records that absolute path and the byte comparison would otherwise differ by the path alone.
# Workload C (the old no-entry default run) additionally gets --legacy_default true: on the branch an empty input_dir stops the default run, by owner decision.
set -uo pipefail
CODE=${CODE:-$(cd "$(dirname "$0")/../.." && pwd)}
DATA=${DATA:-/workspace/Biomarker_pipeline}
P=$DATA
NF=$DATA/env/bin/nextflow
OUT=$(realpath -m "${1:?usage: tests/regression/run_regression.sh <new dir>}")
[[ -e "$OUT" ]] && { echo "ERROR: $OUT exists (never overwritten)"; exit 2; }
mkdir -p "$OUT"
cd "$CODE"
{ echo "start $(date -u '+%F %T') UTC"; echo "CODE $CODE DATA $DATA"; echo "git $(git rev-parse HEAD) branch $(git branch --show-current)"; echo "questions.yaml sha256: code $(sha256sum < $CODE/config/questions.yaml | cut -c1-64) data $(sha256sum < $DATA/config/questions.yaml | cut -c1-64)"; echo "uncommitted: $(git status --porcelain | wc -l)"
  git status --porcelain; echo "lock: $(sha256sum -c --quiet lock/diagnostic_v1/sha256sums_final.txt && echo 7/7 OK || echo MISMATCH)"; } > "$OUT/run_info.txt"

run_A() {
  local R=$OUT/discover; mkdir -p "$R"; cd "$R"
  "$NF" run "$CODE/nextflow/main.nf" -entry DISCOVER -profile conda --question diagnostic --questions_config "$DATA/config/questions.yaml" \
    --counts_matrix "$P/results/04_counts/counts_featurecounts.tsv" \
    --clinical "$P/results_rc1prep_19/discover_regression_20261015/inputs/clinical_19.tsv" \
    --counts_matrix_secondary "$P/results/04_counts/counts_salmon_tximport.tsv" --counts_rds_secondary "$P/results/04_counts/salmon_tximport.rds" \
    --outdir "$R/out" -w "$R/work" -with-trace "$R/trace.tsv" > "$R/run.log" 2>&1
  echo $? > "$R/exit_code.txt"
}
run_B() {
  bash "$CODE/tests/regression/run_diagnostic_branch.sh" "$OUT/diagnostic" > "$OUT/diagnostic.log" 2>&1
  echo $? > "$OUT/diagnostic_exit_code.txt"
}
run_C() {
  local R=$OUT/stub_default; mkdir -p "$R/ref"; cd "$R"
  "$NF" run "$CODE/nextflow/main.nf" -profile conda -stub-run --legacy_default true --base_dir "$DATA" --outdir "$R/out" --ref_dir "$R/ref" -w "$R/work" > "$R/run.log" 2>&1
  echo $? > "$R/exit_code.txt"
}
run_A & run_B & run_C &
wait
echo "end $(date -u '+%F %T') UTC" >> "$OUT/run_info.txt"
echo "A discover exit $(cat "$OUT/discover/exit_code.txt") | B diagnostic exit $(cat "$OUT/diagnostic_exit_code.txt") | C stub exit $(cat "$OUT/stub_default/exit_code.txt")" | tee -a "$OUT/run_info.txt"
