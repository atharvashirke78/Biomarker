#!/usr/bin/env bash
# Regression test for nextflow/bin/nested_full_validation.R (fully nested internal validation). Writes ONLY into a scratch folder, never into the repo.
#   1. 19-sample development matrix, small settings (1 repeat, 3 permutations): runs end to end, the four outputs exist, and TWO runs give byte-identical outputs
#   2. synthetic pure-noise labels: the NESTED AUC is near 0.5, while the same data screened on ALL samples (--screen-outside true) gives a much higher AUC
#   3. usage errors exit non-zero; an analysis problem (too few samples) exits 0 with a 'NOT RUN' row
# usage: tests/regression/test_nested_full_validation.sh [scratch_dir]      env: R_DGE_RSCRIPT, R_CLF_RSCRIPT, COUNTS19, SHEET19, CORES (default 2, max 4)
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
SCRIPT="$REPO/nextflow/bin/nested_full_validation.R"
RS_CLF="${R_CLF_RSCRIPT:-/workspace/Biomarker_pipeline/env/prefix/r_classifier/bin/Rscript}"
COUNTS19="${COUNTS19:-/workspace/Biomarker_pipeline/results/04_counts/counts_featurecounts.tsv}"
SHEET19="${SHEET19:-/workspace/Biomarker_pipeline/results_rc1prep_19/sample_sheet_rc1prep.tsv}"
CORES="${CORES:-2}"; [ "$CORES" -gt 4 ] && CORES=4
SCR="${1:-$(mktemp -d "${TMPDIR:-/tmp}/nested_test.XXXXXX")}"; mkdir -p "$SCR"; SCR="$(cd "$SCR" && pwd)"
case "$SCR" in "$REPO"*|/workspace/Biomarker_pipeline*) echo "scratch folder must be outside the repo and the pipeline folders: $SCR"; exit 2;; esac
fail=0; ok() { echo "PASS: $1"; }; bad() { echo "FAIL: $1"; fail=1; }
common=(--questions "$REPO/config/questions.yaml" --cores "$CORES")
val() { awk -F'\t' -v m="$2" '$1==m {print $3 ? $2 : $2}' "$1" | head -1; }     # value column of a metric
col() { awk -F'\t' -v m="$2" -v c="$3" '$1==m {print $c}' "$1" | head -1; }
echo "scratch: $SCR; cores $CORES"
# ---- 1. 19 development samples, two identical runs
for f in "$COUNTS19" "$SHEET19"; do [ -r "$f" ] || { echo "SKIP: input not readable: $f"; exit 2; }; done
for run in 1 2; do
  s=$(date +%s)
  "$RS_CLF" "$SCRIPT" --counts "$COUNTS19" --sheet "$SHEET19" --rules "$REPO/config/candidate_rules.yaml" "${common[@]}" --repeats 1 --permutations 3 --out "$SCR/dev_run$run" > "$SCR/dev_run$run.log" 2>&1
  rc=$?; e=$(date +%s); echo "19-sample run $run: exit $rc, wall $((e-s)) s"
  [ $rc -eq 0 ] && ok "run $run exit 0" || bad "run $run exit $rc"
done
for f in nested_validation_summary.tsv nested_selection_frequency.tsv nested_permutation_null.tsv NESTED_VALIDATION.md; do
  [ -s "$SCR/dev_run1/$f" ] && ok "output $f exists" || bad "output $f missing"
  cmp -s "$SCR/dev_run1/$f" "$SCR/dev_run2/$f" && ok "$f identical across the two runs" || bad "$f differs between the two runs"
done
grep -q "NOT RUN" "$SCR/dev_run1/nested_validation_summary.tsv" && bad "19-sample run reported NOT RUN" || ok "19-sample run produced an estimate"
grep -q "NESTED, rules 1-2 only; NOT independent validation" "$SCR/dev_run1/nested_validation_summary.tsv" && ok "label present" || bad "label missing"
grep -q "rule 3b" "$SCR/dev_run1/NESTED_VALIDATION.md" && ok "LIMITATION text names rule 3b" || bad "LIMITATION text missing"
echo "19-sample nested AUC: $(col "$SCR/dev_run1/nested_validation_summary.tsv" nested_AUC 2), folds skipped: $(col "$SCR/dev_run1/nested_validation_summary.tsv" folds_skipped 2), permutation p: $(col "$SCR/dev_run1/nested_validation_summary.tsv" permutation_p_value 2)"
# ---- 2. synthetic pure noise: nested vs screening on all samples
SY="$SCR/syn"; mkdir -p "$SY"
cat > "$SY/gen.R" <<'RSCRIPT'
set.seed(11); G <- 1500; n <- 30
batch <- rep(c("A","B"), c(16,14)); group <- c(rep(c("normal","tumour"), c(8,8)), rep(c("normal","tumour"), c(7,7)))
mu <- exp(rnorm(G, 5, 1.5)); bf <- exp(rnorm(G, 0, 0.3)); lib <- exp(rnorm(n, 0, 0.2))
M <- sapply(1:n, function(j) rnbinom(G, mu = mu * lib[j] * (if (batch[j] == "B") bf else 1), size = 4))      # NO group effect at all: the labels are pure noise
rownames(M) <- sprintf("G%04d", 1:G); colnames(M) <- sprintf("S%02d", 1:n)
write.table(data.frame(gene_id = rownames(M), M, check.names = FALSE), "counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample_id = colnames(M), label = colnames(M), group = group, baseline = group == "normal", study = batch, batch = batch, count_source = "featurecounts"), "sheet.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
RSCRIPT
(cd "$SY" && "$RS_CLF" gen.R) || { bad "synthetic generator failed"; exit 1; }
# With a true null DESeq2 + padj < 0.05 finds (almost) nothing, so the demonstration relaxes rule 1 (padj_max 1, |LFC| >= 0.1) in a COPY of the rules file: a permissive screen is what makes selection bias visible.
sed -e 's/padj_max: 0.05/padj_max: 1.0/' -e 's/abs_shrunk_log2fc_min: 1.0/abs_shrunk_log2fc_min: 0.1/' "$REPO/config/candidate_rules.yaml" > "$SY/rules_relaxed.yaml"
for m in false true; do
  s=$(date +%s)
  "$RS_CLF" "$SCRIPT" --counts "$SY/counts.tsv" --sheet "$SY/sheet.tsv" --rules "$SY/rules_relaxed.yaml" "${common[@]}" --repeats 2 --permutations 3 --screen-outside $m --out "$SY/out_$m" > "$SY/out_$m.log" 2>&1
  rc=$?; e=$(date +%s); echo "synthetic (screen-outside=$m): exit $rc, wall $((e-s)) s"; [ $rc -eq 0 ] || bad "synthetic run screen-outside=$m exit $rc"
done
N=$(col "$SY/out_false/nested_validation_summary.tsv" nested_AUC 2); B=$(col "$SY/out_true/nested_validation_summary.tsv" nested_AUC 2)
echo "synthetic pure-noise labels: NESTED AUC = $N ; screened on ALL samples AUC = $B"
awk -v n="$N" 'BEGIN{exit !(n>=0.25 && n<=0.65)}' && ok "nested AUC near 0.5 (0.25-0.65)" || bad "nested AUC not near 0.5: $N"
awk -v n="$N" -v b="$B" 'BEGIN{exit !(b>=0.70 && b-n>=0.15)}' && ok "screen-on-all-samples AUC is much higher (>= 0.70 and >= nested + 0.15)" || bad "no selection-bias gap: nested $N vs all-samples $B"
# ---- 3. usage error and NOT RUN behaviour
"$RS_CLF" "$SCRIPT" --counts "$SY/counts.tsv" > /dev/null 2>&1; [ $? -eq 2 ] && ok "missing options: exit 2" || bad "usage error did not exit 2"
"$RS_CLF" "$SCRIPT" --counts "$SY/counts.tsv" --sheet "$SY/sheet.tsv" --questions "$REPO/config/questions.yaml" --rules "$REPO/config/candidate_rules.yaml" --bogus 1 > /dev/null 2>&1; [ $? -eq 2 ] && ok "unknown option: exit 2" || bad "unknown option did not exit 2"
head -1 "$SY/counts.tsv" | cut -f1-4,10-12 > "$SY/tiny.tsv"; tail -n +2 "$SY/counts.tsv" | cut -f1-4,10-12 >> "$SY/tiny.tsv"
"$RS_CLF" "$SCRIPT" --counts "$SY/tiny.tsv" --sheet "$SY/sheet.tsv" --questions "$REPO/config/questions.yaml" --rules "$REPO/config/candidate_rules.yaml" --out "$SY/out_tiny" > /dev/null 2>&1; rc=$?
[ $rc -eq 0 ] && grep -q "NOT RUN" "$SY/out_tiny/nested_validation_summary.tsv" && ok "6 samples: exit 0 with a NOT RUN row" || bad "tiny cohort: exit $rc / no NOT RUN row"
# nothing may have been written into the repo by this test
echo; [ $fail -eq 0 ] && echo "ALL TESTS PASSED" || echo "SOME TESTS FAILED"; exit $fail
