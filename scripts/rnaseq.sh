#!/usr/bin/env bash
# One-command launcher of the folder entry (plan results_rc1prep_19/phase1_plan_20261016/PHASE1_PLAN.md section 7; owner decisions 9a).
#
#   scripts/rnaseq.sh --mode quant|discover|apply --input_dir DIR [--outdir OUT] [options]
#   scripts/rnaseq.sh --make_metadata --input_dir DIR [--outdir OUT]          (writes metadata.DRAFT.csv into a new run folder, then stops)
#   scripts/rnaseq.sh --mode M --input_dir DIR [--outdir OUT] --run_id <existing run folder name> --resume
#
# DIR holds the FASTQ (<sample>_R1/_R2 or _1/_2, .fastq.gz/.fq.gz; conf/input.config) and metadata.csv (sample_id, group, batch_or_study
# [, patient_id, fastq_R1, fastq_R2, age, ...]). OUT defaults to DIR/result. Every run gets a NEW folder OUT/<mode>_<YYYYMMDD_HHMMSS>/
# (never overwritten); --run_id NAME --resume continues that folder (refused if it is missing or was started in another mode).
# Options:
#   --counts TSV                      counts-level run: featureCounts-style matrix (gene_id + one column per sample); no FASTQ step
#   --counts_secondary TSV --rds_secondary RDS   discover: second route (Salmon/tximport) for rule 3a
#   --qc_table TSV                    apply with --counts: sample QC table (FAIL excluded, WARN kept)
#   --holdout FILE                    discover: hold-out sample IDs (passed to -entry DISCOVER)
#   --cohort development|validation   apply: required, written into every score output (development = optimistic, not validation)
#   --stub                            Nextflow -stub-run for the FASTQ path (structure test; the memory guard still applies)
# The memory guard is never bypassed (no --allow_salmon_only). Nothing is written into lock/ or DIR.
set -uo pipefail
P=$(cd "$(dirname "$0")/.." && pwd)
NF=$P/env/bin/nextflow
PY="$P/env/prefix/pyqc/bin/python3 -I"
RS=$P/env/prefix/r_classifier/bin/Rscript
die() { echo "ERROR: $*" >&2; exit 2; }

MODE= ; IN= ; OUTD= ; RUN_ID= ; RESUME=0 ; COUNTS= ; C2= ; R2= ; QCT= ; HOLD= ; COHORT= ; STUB=0 ; MAKEMETA=0
CMDLINE="$0 $*"
while (( $# )); do
  case "$1" in
    --mode) MODE=${2:-}; shift 2 ;;
    --input_dir) IN=${2:-}; shift 2 ;;
    --outdir) OUTD=${2:-}; shift 2 ;;
    --run_id) RUN_ID=${2:-}; shift 2 ;;
    --resume|-resume) RESUME=1; shift ;;
    --counts) COUNTS=${2:-}; shift 2 ;;
    --counts_secondary) C2=${2:-}; shift 2 ;;
    --rds_secondary) R2=${2:-}; shift 2 ;;
    --qc_table) QCT=${2:-}; shift 2 ;;
    --holdout) HOLD=${2:-}; shift 2 ;;
    --cohort) COHORT=${2:-}; shift 2 ;;
    --stub) STUB=1; shift ;;
    --make_metadata) MAKEMETA=1; MODE=${MODE:-discover}; shift ;;
    -h|--help) sed -n 2,20p "$0"; exit 0 ;;
    *) die "unknown argument: $1 (see --help)" ;;
  esac
done
[[ $MODE =~ ^(quant|discover|apply)$ ]] || die "--mode quant|discover|apply is required"
[[ -n $IN && -d $IN ]] || die "--input_dir must be an existing folder"
IN=$(realpath "$IN")
OUTD=$(realpath -m "${OUTD:-$IN/result}")
abs() { [[ -z $1 ]] && return 0; [[ -e $1 ]] || die "file not found: $1"; realpath "$1"; }
COUNTS=$(abs "$COUNTS") || exit 2; C2=$(abs "$C2") || exit 2; R2=$(abs "$R2") || exit 2
QCT=$(abs "$QCT") || exit 2; HOLD=$(abs "$HOLD") || exit 2
[[ -n $C2 && -z $R2 || -z $C2 && -n $R2 ]] && die "--counts_secondary and --rds_secondary go together"
[[ $MODE == apply && $MAKEMETA == 0 && ! $COHORT =~ ^(development|validation)$ ]] && die "apply needs --cohort development|validation"
[[ $MODE == quant && -n $COUNTS ]] && die "--counts makes no sense with --mode quant"
cd "$P"

# ---- run folder: new, or an existing one with --run_id ... --resume
if (( RESUME )); then
  [[ -n $RUN_ID ]] || die "--resume needs --run_id <existing run folder name>"
  RUN=$OUTD/$RUN_ID
  [[ -d $RUN && -f $RUN/run_info.json ]] || die "run folder not found (nothing to resume): $RUN"
  $PY scripts/run_info.py resume-check --run "$RUN" --mode "$MODE" --input-dir "$IN" || exit 2
  echo "resuming $RUN"
else
  [[ -n $RUN_ID ]] && die "--run_id is only used with --resume (new runs get a new dated folder)"
  RUN=$OUTD/${MODE}_$(date -u +%Y%m%d_%H%M%S)
  [[ -e $RUN ]] && die "$RUN exists (never overwritten); start again in a second"
  mkdir -p "$RUN"
  $PY scripts/run_info.py init --run "$RUN" --mode "$MODE" --input-dir "$IN" --counts "$COUNTS" --cmd "$CMDLINE" || exit 2
  MM=(); (( MAKEMETA )) && MM=(--make-metadata)
  CC=(); [[ -n $COUNTS ]] && CC=(--counts "$COUNTS")
  $PY scripts/input_scan.py --input-dir "$IN" --mode "$MODE" --out-dir "$RUN/input" "${CC[@]}" "${MM[@]}" 2>&1 | tee "$RUN/input_scan.log"
  rc=${PIPESTATUS[0]}
  if (( rc != 0 )); then $PY scripts/run_info.py end --run "$RUN" --status "input_stop" --exit "$rc"; echo "run folder: $RUN"; exit "$rc"; fi
  [[ -f $RUN/input/metadata.used.csv ]] && cp "$RUN/input/metadata.used.csv" "$RUN/metadata.used.csv"
fi
RES=(); (( RESUME )) && RES=(-resume)
STB=(); (( STUB )) && STB=(-stub-run)

# ---- stages; each Nextflow stage has its own launch folder (own .nextflow history, so -resume finds it) and work dir
nf_stage() {  # name, then nextflow args
  local name=$1; shift
  mkdir -p "$RUN/nf_$name"
  echo "[$(date -u '+%F %T') UTC] stage $name"
  (cd "$RUN/nf_$name" && "$NF" run "$P/nextflow/main.nf" -profile conda "$@" -w "$RUN/work_$name" \
      -with-trace "$RUN/nf_$name/trace_$(date -u +%Y%m%d_%H%M%S).tsv" "${RES[@]}" > "$RUN/nf_$name/run_$(date -u +%Y%m%d_%H%M%S).log" 2>&1)
  local rc=$?
  echo "$rc" > "$RUN/nf_$name/exit_code.txt"
  echo "[$(date -u '+%F %T') UTC] stage $name exit $rc"
  return $rc
}
finish() {  # status, exit code
  $PY scripts/run_info.py end --run "$RUN" --status "$1" --exit "$2"
  $PY scripts/make_report.py --run "$RUN" || echo "WARNING: report step failed" >&2
  echo "run folder: $RUN"
  exit "$2"
}

if [[ -z $COUNTS ]]; then
  nf_stage quant -entry QUANT --samplesheet "$RUN/input/samplesheet_quant.csv" --legacy_samplesheet false \
      --outdir "$RUN/pipeline/quant" "${STB[@]}" || finish quant_failed 1
  QDIR=$RUN/pipeline/quant/contract
  CNT=$QDIR/counts_matrix.tsv; C2=$QDIR/counts_salmon_tximport.tsv; R2=$QDIR/salmon_tximport.rds; QCT=$QDIR/qc_table.tsv
  [[ $MODE == quant ]] && finish completed 0
  (( STUB )) && finish "completed (stub: downstream not run on stub counts)" 0
  # samples excluded by QUANT QC (FAIL) are not in the contract counts: keep only clinical rows present in the counts
  $PY scripts/run_info.py filter-clinical --run "$RUN" --counts "$CNT" || finish input_stop 2
  CLIN=$RUN/input/clinical_discover.qcpass.tsv
else
  CNT=$COUNTS; CLIN=$RUN/input/clinical_discover.tsv
fi

if [[ $MODE == discover ]]; then
  SEC=(); [[ -n $C2 ]] && SEC=(--counts_matrix_secondary "$C2" --counts_rds_secondary "$R2")
  HO=(); [[ -n $HOLD ]] && HO=(--holdout "$HOLD")
  nf_stage discover -entry DISCOVER --question diagnostic --counts_matrix "$CNT" --clinical "$CLIN" "${SEC[@]}" "${HO[@]}" \
      --outdir "$RUN/pipeline/discover" || finish discover_failed 1
  finish completed 0
fi

# apply: locked model lock/diagnostic_v1 (hash-checked by the script); evaluation only when every sample has a known label
SC=$RUN/pipeline/apply/scores
if [[ -e $SC ]]; then echo "apply scores exist ($SC): not re-run on resume"; finish completed 0; fi
mkdir -p "$RUN/pipeline/apply"
EV=(); $PY scripts/run_info.py labels-complete --run "$RUN" && EV=(--evaluate)
[[ ${#EV[@]} -gt 0 && $COHORT == validation && -z $QCT ]] && { echo "NOTE: no QC table: evaluation skipped (validation needs --qc_table)"; EV=(); }
QA=(); [[ -n $QCT ]] && QA=(--qc-table "$QCT")
$RS scripts/apply_locked_model.R --counts "$CNT" --sample-sheet "$RUN/input/sample_sheet_apply.tsv" "${QA[@]}" --cohort "$COHORT" "${EV[@]}" \
    --out-dir "$SC" > "$RUN/pipeline/apply/apply.log" 2>&1
rc=$?; cat "$RUN/pipeline/apply/apply.log"
(( rc == 0 )) && finish completed 0 || finish apply_failed 1
