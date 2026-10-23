#!/usr/bin/env bash
# Fully nested internal validation at DEFAULT settings (5 repeats x 5 outer folds, 20 label permutations x 5 repeats) on the 19 development samples.
# Report-only; writes only into the output folder. Meant for an idle pod: it uses most CPUs. About 525 DESeq2 fits (about 7 s each on one core).
# usage: bash scripts/run_nested_default.sh [output folder]      env: NESTED_CORES (default nproc-2, at most 24)
#        run it in tmux:  tmux new -d -s nested 'bash /workspace/atharva_pips/scripts/run_nested_default.sh'
set -uo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
INSTALL=${BIOMARKER_INSTALL_DIR:-/workspace/Biomarker_pipeline}
COUNTS=${COUNTS19:-$INSTALL/results/04_counts/counts_featurecounts.tsv}
SHEET=${SHEET19:-$INSTALL/results_rc1prep_19/sample_sheet_rc1prep.tsv}
OUT=${1:-/workspace/atharva_runs/nested_default_$(date +%Y%m%d_%H%M%S)}
N=$(nproc); CORES=${NESTED_CORES:-$((N > 3 ? N - 2 : 1))}; (( CORES > 24 )) && CORES=24
[[ -e $OUT ]] && { echo "ERROR: $OUT exists (never overwritten)"; exit 2; }
for f in "$COUNTS" "$SHEET" "$INSTALL/env/prefix/r_dge/bin/Rscript" "$INSTALL/env/prefix/r_classifier/bin/Rscript"; do [[ -r $f ]] || { echo "ERROR: not readable: $f"; exit 2; }; done
mkdir -p "$OUT"
export NESTED_MAX_CORES=$CORES OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
{ echo "host $(hostname)  cpus $N  cores used $CORES  load $(cut -d' ' -f1-3 /proc/loadavg)"; echo "start $(date -u '+%F %T') UTC"
  echo "code $(git -C "$HERE" rev-parse --short HEAD 2>/dev/null || echo unknown)  script sha256 $(sha256sum "$HERE/nextflow/bin/nested_full_validation.R" | cut -c1-16)"
  echo "counts $COUNTS"; echo "sheet  $SHEET"; } | tee "$OUT/run_info.txt"
S=$(date +%s)
"$INSTALL/env/prefix/r_classifier/bin/Rscript" "$HERE/nextflow/bin/nested_full_validation.R" --counts "$COUNTS" --sheet "$SHEET" \
  --questions "$HERE/config/questions.yaml" --rules "$HERE/config/candidate_rules.yaml" --thresholds "$HERE/conf/thresholds.config" \
  --cores "$CORES" --out "$OUT" 2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
echo "exit $rc  wall $(( $(date +%s) - S )) s  end $(date -u '+%F %T') UTC" | tee -a "$OUT/run_info.txt"
[[ -s $OUT/nested_validation_summary.tsv ]] && { echo; column -t -s$'\t' "$OUT/nested_validation_summary.tsv" | cut -c1-200; }
exit $rc
