#!/usr/bin/env bash
# One-command diagnostic workflow on the 19 DEVELOPMENT samples from the existing quantification (no FASTQ step).
# Re-runs the analysis stage and applies the LOCKED model; nothing is refitted and no locked file is changed.
#   1. nextflow -entry ANALYSIS: input checks -> DESeq2 (~ batch + group) -> pipeline ORA/GSEA (GO BP, KEGG, Hallmark; 07_pathway).
#      The Reactome analysis of the report is a separate script (scripts/downstream/reactome_analysis.R) and is not re-run here.
#   2. scripts/apply_locked_model.R on the same counts (development data: optimistic, NOT validation)
#   3. scripts/rule4_reproduction.py on the new DESeq2 table (self-consistency only: must reproduce all 1,390 candidates)
#   4. tests/regression/check_19_development.py against the documented expected results
# usage: scripts/run_diagnostic.sh <new output dir>      (the folder must not exist; ~10 min on 8 cores, < 2 GB RAM)
set -euo pipefail
P=/workspace/Biomarker_pipeline
cd "$P"
OUT=${1:?usage: scripts/run_diagnostic.sh <new output dir>}
OUT=$(realpath -m "$OUT")
[[ -e "$OUT" ]] && { echo "ERROR: $OUT exists (never overwritten)"; exit 2; }
COUNTS=$P/results/04_counts/counts_featurecounts.tsv
SHEET=$P/results_rc1prep_19/sample_sheet_rc1prep.tsv
QC=$P/results/04_counts/sample_qc.tsv
mkdir -p "$OUT"
{ echo "start $(date -u '+%F %T') UTC"; echo "git $(git rev-parse --short HEAD) $(git status --porcelain | wc -l) uncommitted"; } > "$OUT/run_info.txt"

echo "[1/4] -entry ANALYSIS"
(cd "$OUT" && "$P/env/bin/nextflow" run "$P/nextflow/main.nf" -entry ANALYSIS -profile conda \
   --counts_matrix "$COUNTS" --sample_sheet "$SHEET" --count_source featurecounts \
   --outdir "$OUT/analysis" -w "$OUT/work" -with-trace "$OUT/trace.tsv" > "$OUT/analysis_run.log" 2>&1)

echo "[2/4] locked model (development)"
env/prefix/r_classifier/bin/Rscript scripts/apply_locked_model.R --counts "$COUNTS" --sample-sheet "$SHEET" --qc-table "$QC" \
  --cohort development --evaluate --out-dir "$OUT/locked_model_development"

echo "[3/4] rule 4 self-consistency"
env/prefix/pyqc/bin/python3 -I scripts/rule4_reproduction.py --results "$OUT/analysis/06_dge/deseq2_results_all.tsv" \
  --summary "$OUT/analysis/06_dge/deseq2_summary.json" --candidates results_rc1prep_19/final_candidates/final_candidates_18of19.tsv \
  --label self-consistency --out-dir "$OUT/rule4_self_consistency"

echo "[4/4] regression checks"
env/prefix/pyqc/bin/python3 -I tests/regression/check_19_development.py --featurecounts-only --de-dir "$OUT/analysis/06_dge" \
  --scores "$OUT/locked_model_development/scores.tsv" --rule4 "$OUT/rule4_self_consistency/rule4_per_gene.tsv" \
  --out "$OUT/regression_summary.tsv"
for f in deseq2_results_all.tsv deseq2_summary.json; do
  cmp -s "$OUT/analysis/06_dge/$f" "results_rc1prep_19/$f" && echo "$f byte-identical to results_rc1prep_19/$f" | tee -a "$OUT/run_info.txt" \
    || { echo "$f DIFFERS from results_rc1prep_19/$f" | tee -a "$OUT/run_info.txt"; exit 1; }
done
echo "end $(date -u '+%F %T') UTC" >> "$OUT/run_info.txt"
echo "done: $OUT"
