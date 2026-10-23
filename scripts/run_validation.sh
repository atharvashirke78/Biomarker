#!/usr/bin/env bash
# Independent validation of the LOCKED diagnostic model (tag lock-diagnostic-v1) on the 15 held-out samples
# (validation_run/samplesheets/validation_15.csv; FASTQs in validation_run/data/{normal,tumor}_sample/).
#
#   scripts/run_validation.sh preflight   checks only (resources, lock hashes, sheet, FASTQ sizes vs ENA); changes nothing
#   scripts/run_validation.sh quant       preflight, then -entry QUANT in the background (nohup); log $OUT/quant_run.log
#   scripts/run_validation.sh resume      resume the QUANT run (same session) after a stop
#   scripts/run_validation.sh downstream  after QUANT finished: -entry ANALYSIS (DESeq2 for rule 4), locked model (applied once), rule 4
#
# OUT defaults to results_rc1prep_15quant (git-ignored). Same settings as the 19-sample run (nextflow.config, conf/thresholds.config);
# no parameter is changed here except paths. The memory guard is never bypassed (no --allow_salmon_only).
set -euo pipefail
P=/workspace/Biomarker_pipeline
cd "$P"
OUT=${OUT:-$P/results_rc1prep_15quant}
SHEET=$P/validation_run/samplesheets/validation_15.csv
NF=$P/env/bin/nextflow
MIN_MEM_BYTES=${MIN_MEM_BYTES:-60000000000}   # what the run needs in practice: STAR ~28.5 GB + Salmon 28 GB x 2 forks (nextflow.config executor budget 58 GB)
MIN_CPUS=${MIN_CPUS:-16}                      # nextflow.config executor cpus = 16
mode=${1:-preflight}

preflight() {
  local ok=1
  local mem; mem=$(cat /sys/fs/cgroup/memory.max 2>/dev/null || cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo unreadable)
  local guard_gib; guard_gib=$(grep -E '^\s*min_mem_star_gb' conf/thresholds.config | grep -oE '[0-9]+' | head -1)
  echo "cgroup memory limit: $mem bytes | pipeline guard min_mem_star_gb: ${guard_gib} GiB | required for this run: >= $MIN_MEM_BYTES bytes"
  if ! [[ "$mem" =~ ^[0-9]+$ ]] || (( mem < MIN_MEM_BYTES )); then echo "BLOCKED: memory below requirement"; ok=0; fi
  echo "cpus: $(nproc) | required: >= $MIN_CPUS"
  if (( $(nproc) < MIN_CPUS )); then echo "BLOCKED: fewer cores than required"; ok=0; fi
  echo "free disk on /workspace: $(df -h /workspace | awk 'NR==2{print $4}') (need about 4 x 96 GB FASTQ = ~400 GB for BAM + Salmon working files)"
  if sha256sum -c --quiet lock/diagnostic_v1/sha256sums_final.txt; then echo "lock hashes: OK"; else echo "BLOCKED: lock hashes do not match"; ok=0; fi
  env/prefix/pyqc/bin/python3 -I - "$SHEET" <<'EOF' || ok=0
import sys, os, pandas as pd
s = pd.read_csv(sys.argv[1])
meta = pd.concat([pd.read_csv(f'results/10_independent_validation/{st}/metadata/ena_filereport.tsv', sep='\t') for st in ('SRP344867', 'SRP245232')]).set_index('run_accession')
bad = []
for _, r in s.iterrows():
    d = 'validation_run/data/' + ('normal_sample' if r.group == 'normal' else 'tumor_sample')
    sizes = [int(x) for x in str(meta.loc[r.run, 'fastq_bytes']).split(';')]
    for i, e in enumerate(('R1', 'R2')):
        p = f'{d}/{r.run}_{e}.fastq.gz'
        if not os.path.exists(p) or os.path.getsize(p) != sizes[i]:
            bad.append(p)
dev = set(pd.read_csv('samplesheets/diagnostic_20.csv').run)
print(f'sample sheet: {len(s)} rows, unique {s["sample"].is_unique}, overlap with development runs {len(dev & set(s.run))}; FASTQ files present with ENA size {2 * len(s) - len(bad)}/{2 * len(s)}')
sys.exit(1 if bad or not s['sample'].is_unique or dev & set(s.run) else 0)
EOF
  if (( ok )); then echo "PREFLIGHT: PASS"; return 0; else echo "PREFLIGHT: FAIL (nothing was started)"; return 1; fi
}

case "$mode" in
  preflight) preflight ;;
  quant)
    preflight
    [[ -e "$OUT/quant" ]] && { echo "ERROR: $OUT/quant exists; use 'resume' or another OUT"; exit 2; }
    mkdir -p "$OUT"; cd "$OUT"
    nohup "$NF" run "$P/nextflow/main.nf" -entry QUANT -profile conda \
      --base_dir "$P/validation_run" --samplesheet "$SHEET" --outdir "$OUT/quant" -w "$OUT/work" \
      -with-trace "$OUT/quant_trace.tsv" -with-report "$OUT/quant_report.html" > "$OUT/quant_run.log" 2>&1 &
    echo "QUANT started (pid $!). Watch: tail -f $OUT/quant_run.log   Resume: OUT=$OUT scripts/run_validation.sh resume" ;;
  resume)
    preflight
    cd "$OUT"
    nohup "$NF" run "$P/nextflow/main.nf" -entry QUANT -profile conda -resume \
      --base_dir "$P/validation_run" --samplesheet "$SHEET" --outdir "$OUT/quant" -w "$OUT/work" \
      -with-trace "$OUT/quant_trace_resume_$(date +%Y%m%d_%H%M%S).tsv" >> "$OUT/quant_run.log" 2>&1 &
    echo "QUANT resumed (pid $!). Watch: tail -f $OUT/quant_run.log" ;;
  downstream)
    C="$OUT/quant/contract"
    for f in counts_matrix.tsv sample_sheet.tsv qc_table.tsv; do [[ -s "$C/$f" ]] || { echo "ERROR: $C/$f missing: QUANT has not finished"; exit 2; }; done
    [[ -e "$OUT/analysis" || -e "$OUT/locked_model_validation" || -e "$OUT/rule4" ]] && { echo "ERROR: downstream outputs exist; the locked model is applied once only"; exit 2; }
    (cd "$OUT" && "$NF" run "$P/nextflow/main.nf" -entry ANALYSIS -profile conda --base_dir "$P" \
       --counts_matrix "$C/counts_matrix.tsv" --sample_sheet "$C/sample_sheet.tsv" --qc_table "$C/qc_table.tsv" \
       --quant_manifest "$C/manifest_quant.json" --outdir "$OUT/analysis" -w "$OUT/work_analysis" > "$OUT/analysis_run.log" 2>&1)
    env/prefix/r_classifier/bin/Rscript scripts/apply_locked_model.R --counts "$C/counts_matrix.tsv" --sample-sheet "$C/sample_sheet.tsv" \
      --qc-table "$C/qc_table.tsv" --cohort validation --evaluate --out-dir "$OUT/locked_model_validation"
    env/prefix/pyqc/bin/python3 -I scripts/rule4_reproduction.py --results "$OUT/analysis/06_dge/deseq2_results_all.tsv" \
      --summary "$OUT/analysis/06_dge/deseq2_summary.json" --candidates results_rc1prep_19/final_candidates/final_candidates_18of19.tsv \
      --label independent --out-dir "$OUT/rule4"
    echo "downstream done: $OUT/locked_model_validation, $OUT/rule4, $OUT/analysis" ;;
  *) echo "usage: $0 preflight|quant|resume|downstream"; exit 2 ;;
esac
