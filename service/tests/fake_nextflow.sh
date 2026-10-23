#!/usr/bin/env bash
# Stand-in for `nextflow run main.nf ...` used only by the service tests: records its arguments and writes a minimal output folder.
# Behaviour switches (files in the input folder): FAIL -> input-check style failure (exit 1); NOREPORT -> exit 0 without a report.
while [[ $# -gt 0 ]]; do case "$1" in --input_dir) IN="$2"; shift 2;; --output_dir) OUT="$2"; shift 2;; *) ARGS+=("$1"); shift;; esac; done
mkdir -p "$OUT/logs" "$OUT/10_report/report" "$OUT/10_report/plots" "$OUT/work/ab"
echo "input_dir=$IN output_dir=$OUT extra=${ARGS[*]}" > "$OUT/argv.txt"
echo "secret" > "$OUT/work/ab/secret.txt"
# like the real thing: the trace lists finished tasks only; submitted tasks appear in the console
printf 'task_id\thash\tnative_id\tname\tstatus\texit\n1\taa/bb\t1\tRUN:RUN_PIPELINE:QUANT_STAGE:FASTP (S1)\tCOMPLETED\t0\n' > "$OUT/logs/trace_20260101_000000.tsv"
echo '[aa/bb] Submitted process > RUN:RUN_PIPELINE:QUANT_STAGE:FASTP (S1)'
echo '[cc/dd] Submitted process > RUN:RUN_PIPELINE:QUANT_STAGE:SALMON_QUANT (S1)'
sleep "${FAKE_SLEEP:-1}"
if [[ -e "$IN/FAIL" ]]; then echo "INPUT CHECK FAILED (nothing was run). boom" ; exit 1; fi
[[ -e "$IN/NOREPORT" ]] && exit 0
mkdir -p "$OUT/04_counts" "$OUT/09_discover/02_candidates" "$OUT/09_discover/03_model"
printf 'srr\tlabel\tstudy\tgroup\tstar_unique_pct\tfc_assigned_pct\tsalmon_mapped_pct\tfastp_q30_pct\tfastp_pass_pct\tstatus\treason\nS1\tS1\tB1\ttumour\t85\t60\t90\t95\t98\tPASS\t\nS2\tS2\tB1\tnormal\t60\t30\t90\t95\t98\tWARN\tfc assigned 30%% < 60\n' > "$OUT/04_counts/sample_qc.tsv"
printf 'step\tall_genes\tprotein_coding\nrule 1 (single route)\t100\t90\n+ rule 3b\t10\t9\n' > "$OUT/09_discover/02_candidates/funnel_counts.tsv"
printf 'gene_id\tsymbol\tall_rules_pass\tlog2FC_shrunk_fc\tpadj_fc\tloo_n_pass\nENSG1\tGENEA\tTRUE\t2.1\t0.001\t6\nENSG2\tGENEB\tFALSE\t1.2\t0.04\t3\n' > "$OUT/09_discover/02_candidates/candidate_table.tsv"
printf 'gene_id\tsymbol\tset\tlasso_coef_lambda1se\nENSG1\tGENEA\tup\t0.5\n' > "$OUT/09_discover/03_model/model_genes.tsv"
printf 'model: balanced\ncutoff = 1.25   (median of 5 inner-fold cutoffs)\n' > "$OUT/09_discover/03_model/model_spec.txt"
printf 'n\tAUC\tci_low\tci_high\tTP\tFN\tTN\tFP\tlabel\n6\t1\t1\t1\t3\t0\t3\t0\tdevelopment data: OPTIMISTIC\n' > "$OUT/09_discover/03_model/development_performance.tsv"
echo "<html>report <img src=\"../plots/p.png\"></html>" > "$OUT/10_report/report/report.html"; echo png > "$OUT/10_report/plots/p.png"; echo png > "$OUT/10_report/plots/deseq2_volcano.png"
exit 0
