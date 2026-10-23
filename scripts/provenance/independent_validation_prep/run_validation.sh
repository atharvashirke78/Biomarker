#!/usr/bin/env bash
# Launch the unchanged pipeline (nextflow/main.nf) on the 15 validation samples from validation_run/ with its own work dir and outdir.
# No -resume of the main session; main .nextflow/history is not touched (launch dir = validation_run).
# Runs from a plain shell: project env bin paths are exported below. Aborts (exit 3) on a pod that is too small, so the
# pipeline can never silently fall back to another method (main.nf turns STAR off by itself when the cgroup limit is < 32 GB).
B=/workspace/Biomarker_pipeline; V=$B/validation_run
export PATH=$B/env/bin:$B/env/mamba/envs/nf/bin:$B/env/prefix/r_dge/bin:$B/env/prefix/star/bin:$B/env/prefix/salmon/bin:$B/env/prefix/fastp/bin:$B/env/prefix/samtools/bin:$B/env/prefix/subread/bin:$PATH
MIN_MEM=60000000000; MIN_CPU=16
TS=$(date +%Y%m%d_%H%M%S); mkdir -p $B/logs $B/results/10_independent_validation/all15
LOG=$B/logs/validation_run_$TS.log
say() { echo "$@" | tee -a "$LOG"; }

# ---- preflight
MEM=$(cat /sys/fs/cgroup/memory.max 2>&1); CPU=$(nproc)
say "[preflight] $(date '+%F %T')  memory.max=$MEM  nproc=$CPU  (required: memory.max >= $MIN_MEM, nproc >= $MIN_CPU)"
case "$MEM" in
  max) say "[preflight] memory.max is 'max' (no cgroup limit): accepted, but main.nf then cannot see a limit; check the pod by hand" ;;
  ''|*[!0-9]*) say "[preflight] ABORT: cannot read a numeric memory.max ('$MEM')"; exit 3 ;;
  *) if [ "$MEM" -lt "$MIN_MEM" ]; then
       say "[preflight] ABORT: memory.max $MEM bytes < $MIN_MEM (60 GB). STAR (33 GB) and Salmon (28 GB x2) do not fit; main.nf would disable STAR and change the method. Use a 64 GB pod."; exit 3
     fi ;;
esac
if [ "$CPU" -lt "$MIN_CPU" ]; then say "[preflight] ABORT: nproc $CPU < $MIN_CPU (config assumes 16 cores). Use a 16-core pod."; exit 3; fi
for t in nextflow STAR salmon fastp R; do command -v $t > /dev/null || { say "[preflight] ABORT: $t not found on PATH"; exit 3; }; done
say "[preflight] OK"
say "[versions] $(nextflow -version 2>&1 | grep -i 'version')"
say "[versions] STAR $(STAR --version 2>&1 | head -1)"
say "[versions] $(salmon --version 2>&1 | head -1)"
say "[versions] fastp $(fastp --version 2>&1 | head -1)"
say "[versions] $(R --version 2>&1 | head -1)"

echo "$TS" > $B/logs/validation_latest_ts.txt
cd $V || exit 1
nextflow run $B/nextflow/main.nf -profile conda \
  -w /workspace/validation_work \
  --base_dir $V \
  --samplesheet $V/samplesheets/validation_15.csv \
  --outdir $B/results/10_independent_validation/all15 \
  -with-report   $B/logs/validation_report_$TS.html \
  -with-timeline $B/logs/validation_timeline_$TS.html \
  -with-trace    $B/logs/validation_trace_$TS.tsv \
  "$@" 2>&1 | tee -a "$LOG"
