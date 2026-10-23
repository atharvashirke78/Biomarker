#!/usr/bin/env bash
# Structure test of the one-workflow entry (nextflow/main.nf --input_dir ...) with -stub-run: tiny fake FASTQ, fake reference, no tool is run,
# no real data, no STAR. Needs the install (env/bin/nextflow, env/prefix/pyqc) and a pod that passes the memory guard.
#
# usage: tests/integration/one_workflow_stub.sh --work <NEW folder> [--install-dir /workspace/Biomarker_pipeline]
# Every case prints PASS or FAIL; exit 0 only if all pass. Writes only into --work.
set -uo pipefail
HERE=$(cd "$(dirname "$0")/../.." && pwd)
WORK= ; INSTALL=${BIOMARKER_INSTALL_DIR:-/workspace/Biomarker_pipeline}
while (( $# )); do case "$1" in --work) WORK=${2:-}; shift 2 ;; --install-dir) INSTALL=${2:-}; shift 2 ;; *) echo "unknown argument $1"; exit 2 ;; esac; done
[[ -n $WORK && ! -e $WORK ]] || { echo "usage: $0 --work <new folder> [--install-dir DIR]"; exit 2; }
mkdir -p "$WORK"; WORK=$(cd "$WORK" && pwd)
NF=$INSTALL/env/bin/nextflow; MAIN=$HERE/nextflow/main.nf
[[ -x $NF ]] || { echo "nextflow not found: $NF"; exit 2; }
strip() { sed 's/\x1b\[[0-9;]*m//g'; }
PASS=0; FAIL=0
check() { if eval "$2"; then echo "PASS  $1"; PASS=$((PASS+1)); else echo "FAIL  $1   [$2]"; FAIL=$((FAIL+1)); fi; }

mkfq() { local d=$1; shift; for s in "$@"; do for m in R1 R2; do echo x | gzip > "$d/${s}_$m.fastq.gz"; done; done; }
META_OK=$'sample_id,group,patient_id\na1,normal,P1\na2,normal,P2\na3,normal,P3\nt1,tumor,P1\nt2,tumor,P2\nt3,tumor,P3\n'
newin() { local d=$WORK/in_$1; mkdir -p "$d"; mkfq "$d" ${2:-a1 a2 a3 t1 t2 t3}; [[ ${3-x} != NONE ]] && printf '%s' "${3:-$META_OK}" > "$d/metadata.csv"; echo "$d"; }
# run <launch folder> <log> args...   (stub, fake reference)
run() { local L=$1 log=$2; shift 2; mkdir -p "$L"; (cd "$L" && timeout 900 "$NF" run "$MAIN" -stub-run --ref_dir "$WORK/ref_stub" "$@" 2>&1 | strip > "$log"); return "${PIPESTATUS[0]}"; }
msg() { grep -q -- "$2" "$1"; }

echo "== work folder: $WORK ; install: $INSTALL ; code: $HERE"

# ---- 0. prime the fake reference (a stub run on an empty --ref_dir builds it; later runs find it, as in the real install, so -resume can reuse tasks)
PRIME=$(newin prime); run "$WORK/l0" "$WORK/0.log" --input_dir "$PRIME" --output_dir "$WORK/out_prime"; rc=$?
check "prime: fake reference built" "[[ $rc == 0 && -e $WORK/ref_stub ]]"

# ---- 1. discover, default output folder <input_dir>/Result, then resume
IN=$(newin ok); run "$WORK/l1" "$WORK/1.log" --input_dir "$IN"; rc=$?
R=$IN/Result
check "discover: exit 0" "[[ $rc == 0 ]]"
check "discover: default output folder is <input_dir>/Result" "[[ -d $R/00_input ]]"
for d in 01_fastp 02_alignment 03_quant_salmon 04_counts 05_cohort_qc 06_dge 07_pathway 08_multiqc 09_discover/02_candidates 09_discover/03_model 10_report/report contract logs work; do
  check "discover: $d exists" "[[ -e $R/$d ]]"; done
for f in run_info.json params.json tool_versions.tsv RESULTS_INDEX_RUN.md; do check "discover: $f exists" "[[ -s $R/$f || $f == tool_versions.tsv && -e $R/$f ]]"; done
check "discover: lock 7/7 at start and end" "python3 -c \"import json,sys; a=json.load(open('$R/run_info.json')); sys.exit(0 if a['lock_diagnostic_v1']=='7/7 OK' and a['attempts'][-1]['lock_check_end'].startswith('OK (7 of 7)') else 1)\""
check "discover: Result/ is not scanned as input" "! grep -q 'Result' $R/00_input/problem_files.tsv"
run "$WORK/l1" "$WORK/1b.log" --input_dir "$IN" -resume; rc=$?
tr_=$(ls "$R"/logs/trace_*.tsv | tail -1)
cached=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $h["status"]=="CACHED"{n++} END{print n+0}' "$tr_")
check "resume: exit 0" "[[ $rc == 0 ]]"
check "resume: most tasks cached ($cached)" "[[ $cached -ge 60 ]]"
check "resume: two attempts recorded" "python3 -c \"import json,sys; a=json.load(open('$R/run_info.json'))['attempts']; sys.exit(0 if [x['kind'] for x in a]==['new','resume'] and all(x['status']=='completed' for x in a) else 1)\""

# ---- 2. restart rules
run "$WORK/l1" "$WORK/2a.log" --input_dir "$IN"; rc=$?
check "second start without -resume stops" "[[ $rc != 0 ]] && msg $WORK/2a.log 'already holds a run'"
printf '%s' "${META_OK/P3/P9}" > "$IN/metadata.csv"
run "$WORK/l1" "$WORK/2b.log" --input_dir "$IN" -resume; rc=$?
check "changed input on -resume stops" "[[ $rc != 0 ]] && msg $WORK/2b.log 'input folder changed'"

# ---- 3. apply mode
IN3=$(newin apply); run "$WORK/l3" "$WORK/3.log" --input_dir "$IN3" --output_dir "$WORK/out_apply" --mode apply --cohort development; rc=$?
check "apply: exit 0" "[[ $rc == 0 ]]"
check "apply: 09_locked_model/scores exists" "[[ -d $WORK/out_apply/09_locked_model/scores ]]"
check "apply: no discover folder" "[[ ! -e $WORK/out_apply/09_discover ]]"
IN3b=$(newin apply_nocohort); run "$WORK/l3b" "$WORK/3b.log" --input_dir "$IN3b" --output_dir "$WORK/out_apply_nc" --mode apply; rc=$?
check "apply without --cohort stops" "[[ $rc != 0 ]] && msg $WORK/3b.log 'needs --cohort'"

# ---- 4. stops before any compute (message and no task run)
stopcase() { # name, folder, message fragment, extra args...
  local n=$1 d=$2 frag=$3; shift 3; run "$WORK/ls_$n" "$WORK/s_$n.log" --input_dir "$d" --output_dir "$WORK/out_$n" "$@"; local rc=$?
  check "stop: $n" "[[ $rc != 0 ]] && msg $WORK/s_$n.log '$frag' && [[ ! -e $WORK/out_$n/01_fastp ]]"; }
stopcase no_metadata "$(newin nometa '' NONE)" 'metadata.csv not found'
D=$(newin unpaired); rm "$D/t3_R2.fastq.gz"; stopcase unpaired "$D" 'problem file'
stopcase fastq_without_row "$(newin norow 'a1 a2 a3 t1 t2 t3' "$(head -6 <<<"$META_OK")"$'\n')" 'without a metadata row'
stopcase row_without_fastq "$(newin norowfq 'a1 a2 a3 t1 t2')" 'without FASTQ'
stopcase no_group_column "$(newin nogroup 'a1 a2 a3 t1 t2 t3' $'sample_id,condition,patient_id\na1,normal,P1\n')" 'lacks required column'
stopcase unknown_label "$(newin badlabel 'a1 a2 a3 t1 t2 t3' "${META_OK/t3,tumor/t3,mutant}")" 'unrecognised group label'
stopcase too_few_per_group "$(newin few 'a1 a2 t1 t2 t3' $'sample_id,group,patient_id\na1,normal,P1\na2,normal,P2\nt1,tumor,P1\nt2,tumor,P2\nt3,tumor,P3\n')" 'too few samples in group'

# ---- 5. patient_id is optional (discover): samples are treated as independent, noted in the input summary
NP=$(newin nopatient 'a1 a2 a3 t1 t2 t3' $'sample_id,group\na1,normal\na2,normal\na3,normal\nt1,tumor\nt2,tumor\nt3,tumor\n')
run "$WORK/l5" "$WORK/5.log" --input_dir "$NP"; rc=$?
check "no patient_id: discover runs (exit 0)" "[[ $rc == 0 ]]"
check "no patient_id: noted in input_summary.json" "grep -q 'treated as an independent patient' $NP/Result/00_input/input_summary.json"

# ---- 6. QC stop rule (check_qc_min.py): excluded samples listed, stop with a clear message when too few are left
QCD=$WORK/qc_unit; mkdir -p "$QCD"; HDR=$'srr\tlabel\tstudy\tgroup\tstar_unique_pct\tfc_assigned_pct\tsalmon_mapped_pct\tfastp_q30_pct\tfastp_pass_pct\tstatus\treason'
row() { printf '%s\t%s\ts\t%s\t1\t1\t1\t1\t1\t%s\t%s\n' "$1" "$1" "$2" "$3" "$4"; }
{ echo "$HDR"; row A normal PASS ''; row B tumour FAIL 'Salmon mapped 46% < FAIL limit 50%'; row C normal WARN low; row D tumour PASS ''; } > "$QCD/q1.tsv"
{ echo "$HDR"; row A normal PASS ''; row B normal PASS ''; row C normal PASS ''; row D tumour PASS ''; row E tumour PASS ''; row F tumour PASS ''; } > "$QCD/q2.tsv"
qc() { (cd "$QCD" && "$INSTALL/env/prefix/pyqc/bin/python3" -I "$HERE/nextflow/bin/check_qc_min.py" "$@" 2> err.txt); }
qc --qc q1.tsv --min-total 3; check "qc: 3 of 4 kept (apply) passes and lists the FAIL sample" "[[ $? == 0 ]] && grep -q 'B' $QCD/qc_excluded.tsv"
qc --qc q1.tsv --min-total 3 --min-per-group 3; rc=$?
check "qc: discover stops, names group and excluded sample" "[[ $rc == 2 ]] && grep -q \"group 'tumour' has 1\" $QCD/err.txt && grep -q 'B (tumour)' $QCD/err.txt"
qc --qc q2.tsv --min-total 3 --min-per-group 3; check "qc: 3 + 3 samples pass" "[[ $? == 0 ]]"
head -2 "$QCD/q1.tsv" > "$QCD/q3.tsv"; qc --qc q3.tsv --min-total 3; rc=$?
check "qc: one sample stops with the minimum named" "[[ $rc == 2 ]] && grep -q 'at least 3' $QCD/err.txt"

# ---- 7. no input folder / mistyped input folder: plain message, nothing created
run "$WORK/l7a" "$WORK/7a.log"; rc=$?
check "no input_dir (default workflow) stops with a plain message" "[[ $rc != 0 ]] && msg $WORK/7a.log 'NO INPUT FOLDER SET'"
run "$WORK/l7b" "$WORK/7b.log" -entry RUN; rc=$?
check "no input_dir (-entry RUN) stops at the startup checks" "[[ $rc != 0 ]] && msg $WORK/7b.log 'STARTUP CHECKS FAILED' && msg $WORK/7b.log 'NO INPUT FOLDER SET'"
run "$WORK/l7c" "$WORK/7c.log" --input_dir "$WORK/typo_folder"; rc=$?
check "mistyped input_dir: 'does not exist' and no folder created" "[[ $rc != 0 ]] && msg $WORK/7c.log 'does not exist or is not a folder' && [[ ! -e $WORK/typo_folder ]]"

# ---- 8. startup checks against an empty install (--check_only runs the software and reference checks even under -stub-run)
mkdir -p "$WORK/empty_install"; NI=$(newin checkonly)
(cd "$WORK" && BIOMARKER_INSTALL_DIR="$WORK/empty_install" timeout 300 "$NF" run "$MAIN" -stub-run --check_only true --input_dir "$NI" 2>&1 | strip > "$WORK/8a.log"); rc=${PIPESTATUS[0]}
check "check_only on an empty install lists the missing parts" "[[ $rc != 0 ]] && msg $WORK/8a.log 'micromamba' && msg $WORK/8a.log 'conda environments' && msg $WORK/8a.log 'reference' && msg $WORK/8a.log 'STAR index'"
(cd "$WORK" && timeout 300 "$NF" run "$MAIN" -stub-run --check_only true --input_dir "$NI" 2>&1 | strip > "$WORK/8b.log"); rc=${PIPESTATUS[0]}
check "check_only on the real install: exit 0, 0 missing, no results folder work done" "[[ $rc == 0 ]] && msg $WORK/8b.log '0 missing' && [[ ! -e $NI/Result/10_report ]]"

# ---- 9. --stop_after: stops, writes a QC report, continues with -resume
SA=$(newin stopafter); run "$WORK/l9" "$WORK/9a.log" --input_dir "$SA" --stop_after quant; rc=$?
check "stop_after quant (discover): exit 0, report written, no DE results, no discovery" "[[ $rc == 0 && -d $SA/Result/10_report/report && ! -e $SA/Result/06_dge/deseq2_results_all.tsv && ! -e $SA/Result/09_discover ]]"
check "stop_after quant: recorded in run_info.json" "python3 -c \"import json,sys; sys.exit(0 if json.load(open('$SA/Result/run_info.json'))['attempts'][-1].get('stop_after')=='quant' else 1)\""
run "$WORK/l9" "$WORK/9b.log" --input_dir "$SA" -resume; rc=$?
check "continue after stop_after quant (same command + -resume): runs the rest" "[[ $rc == 0 && -d $SA/Result/06_dge && -d $SA/Result/09_discover/03_model ]]"
SB=$(newin stopapply); run "$WORK/l9c" "$WORK/9c.log" --input_dir "$SB" --mode apply --cohort development --stop_after quant; rc=$?
check "apply + stop_after quant: no locked-model scoring" "[[ $rc == 0 && ! -e $SB/Result/09_locked_model && -d $SB/Result/10_report ]]"
run "$WORK/l9c" "$WORK/9d.log" --input_dir "$SB" --mode apply --cohort development -resume; rc=$?
check "apply: scoring runs after the continue" "[[ $rc == 0 && -d $SB/Result/09_locked_model/scores ]]"
SC=$(newin stopbad); run "$WORK/l9e" "$WORK/9e.log" --input_dir "$SC" --mode apply --cohort development --stop_after analysis; rc=$?
check "apply + stop_after analysis stops with the reason" "[[ $rc != 0 ]] && msg $WORK/9e.log 'applies to --mode discover only'"
run "$WORK/l9f" "$WORK/9f.log" --input_dir "$SC" --stop_after nonsense; rc=$?
check "invalid --stop_after value stops" "[[ $rc != 0 ]] && msg $WORK/9f.log 'allowed values are quant, analysis, discover'"

# ---- 10. the real report script on a staged layout: notices, plot links, name clashes (stub REPORT does not run it)
RL=$WORK/report_layout; R1=$IN/Result; mkdir -p $RL/input $RL/pipeline/quant/04_counts $RL/pipeline/quant/05_cohort_qc $RL/pipeline/analysis/06_dge $RL/pipeline/analysis/07_pathway $RL/pipeline/discover/01_de $RL/pipeline/discover/02_candidates $RL/pipeline/discover/03_model
python3 -c "import json; d=json.load(open('$R1/run_info.json')); d['attempts'][-1]['status']='completed'; d['stop_after']=None; json.dump(d,open('$RL/run_info.json','w'))"
cp $R1/00_input/{scan_files.tsv,problem_files.tsv,metadata.used.csv,metadata.resolved.tsv,group_mapping.tsv,input_summary.json} $RL/input/
printf 'sample\tgroup\tstatus\treason\nX1\ttumour\tFAIL\tSalmon mapped 46%% < FAIL limit 50%%\n' > $RL/input/qc_excluded.tsv
cp $R1/contract/qc_table.tsv $RL/pipeline/quant/04_counts/sample_qc.tsv
for d in 01_de 02_candidates 03_model; do cp $R1/09_discover/$d/* $RL/pipeline/discover/$d/; done
for n in pca_group_study sample_correlation_heatmap library_size_barplot; do printf '\x89PNG' > $RL/pipeline/quant/05_cohort_qc/$n.png; done
for n in deseq2_volcano deseq2_MA deseq2_top50_heatmap; do printf '\x89PNG' > $RL/pipeline/analysis/06_dge/$n.png; printf '\x89PNG' > $RL/pipeline/discover/01_de/$n.png; done
for n in ora_up_GO_BP_dotplot gsea_GO_BP_top; do printf '\x89PNG' > $RL/pipeline/analysis/07_pathway/$n.png; done
mkdir -p $RL/pipeline/quant; printf '{"parameters": {"gencode_release": "v50"}, "references": {"gtf": {"path": "x.gtf.gz", "size": 1, "md5": "abc123"}}}' > $RL/pipeline/quant/manifest_quant.json
printf '{"design": "~ group"}' > $RL/pipeline/discover/01_de/deseq2_summary.json
(cd $RL && "$INSTALL/env/prefix/pyqc/bin/python3" -I "$HERE/scripts/make_report.py" --run . > "$WORK/10.log" 2>&1); rc=$?
check "report: runs on the staged layout" "[[ $rc == 0 && -s $RL/report/report.md && -s $RL/report/report.html ]]"
check "report: mechanics-test notice (smallest group below 10)" "grep -q 'MECHANICS TEST ONLY' $RL/report/report.md"
check "report: QC-excluded samples listed with the reason" "grep -q 'X1 (tumour): Salmon mapped' $RL/report/report.md"
check "report: links volcano, MA, heatmap, pathway and cohort-QC plots" "for k in volcano _MA top50_heatmap dotplot gsea_GO_BP_top pca_group_study sample_correlation_heatmap library_size_barplot; do grep -q \"!\\[.*\\](../plots/.*\$k\" $RL/report/report.md || exit 1; done"
check "report: same-named plots from two stages are both kept" "[[ -e $RL/plots/deseq2_volcano.png && -e $RL/plots/01_de_deseq2_volcano.png ]]"
python3 -c "import json; d=json.load(open('$RL/run_info.json')); d['stop_after']='quant'; json.dump(d,open('$RL/run_info.json','w'))"
sleep 1; (cd $RL && "$INSTALL/env/prefix/pyqc/bin/python3" -I "$HERE/scripts/make_report.py" --run . > "$WORK/10b.log" 2>&1)
check "report: says where the run was stopped (--stop_after)" "grep -q 'stopped after stage .quant.' $RL/report_*/report.md || grep -q 'stopped after stage .quant.' $RL/report/report.md"
check "report: HTML embeds the plots" "[[ \$(grep -c '<img' $RL/report/report.html) -ge 10 ]]"

# ---- 11. report banner follows the cohort (development = NOT the independent validation; validation = no "optimistic" text)
for c in development validation; do
  CL=$WORK/report_cohort_$c; mkdir -p $CL/input $CL/pipeline/apply/scores
  python3 -c "import json; d=json.load(open('$R1/run_info.json')); d['attempts'][-1]['status']='completed'; json.dump(d,open('$CL/run_info.json','w'))"
  cp $R1/00_input/input_summary.json $R1/00_input/group_mapping.tsv $CL/input/
  printf 'sample\tscore\tcall_tumour\tqc_warn\tcohort\ttruth\ttruth_positive\tstudy\nS1\t1.0\t1\tFALSE\t%s\ttumour\t1\tA\n' "$c" > $CL/pipeline/apply/scores/scores.tsv
  (cd $CL && "$INSTALL/env/prefix/pyqc/bin/python3" -I "$HERE/scripts/make_report.py" --run . > "$WORK/11_$c.log" 2>&1)
done
check "report banner (development): says NOT the independent validation" "grep -q 'COHORT: DEVELOPMENT.*NOT the independent validation' $WORK/report_cohort_development/report/report.md"
check "report banner (validation): says validation and has no 'OPTIMISTIC'" "grep -q 'COHORT: VALIDATION' $WORK/report_cohort_validation/report/report.md && ! grep -q 'OPTIMISTIC' $WORK/report_cohort_validation/report/report.md"

# ---- 12. fixes 3, 4, 5 (each with its own PASS line) and 6, 7, 2
# fix 3: group confounded with batch (all normal in batch A, all tumour in batch B) stops at the START, with a plain message, before any compute
stopcase fix3_confounded_batch "$(newin confounded 'a1 a2 a3 t1 t2 t3' $'sample_id,group,batch_or_study,patient_id\na1,normal,A,P1\na2,normal,A,P2\na3,normal,A,P3\nt1,tumor,B,P4\nt2,tumor,B,P5\nt3,tumor,B,P6\n')" 'group is confounded with batch'
# control: both batches contain both groups -> the scan accepts it (fully crossed except the LOO refits, which are skipped and reported later)
CTRL=$(newin crossed 'a1 a2 a3 t1 t2 t3' $'sample_id,group,batch_or_study,patient_id\na1,normal,A,P1\na2,normal,A,P2\na3,normal,B,P3\nt1,tumor,A,P4\nt2,tumor,B,P5\nt3,tumor,B,P6\n')
run "$WORK/l12c" "$WORK/12c.log" --input_dir "$CTRL"; rc=$?
check "fix 3 (control): a batch design with both groups in every batch is accepted" "[[ $rc == 0 ]]"
# fix 4: spaces, commas and quotes in sample IDs stop at the start with a plain message
for nm in space comma quote; do case $nm in space) id='a 1';; comma) id='a,1';; quote) id='a"1';; esac
  d=$WORK/in_badid_$nm; mkdir -p "$d"; mkfq "$d" "$id" a2 a3 t1 t2 t3; printf 'sample_id,group\n"%s",normal\na2,normal\na3,normal\nt1,tumor\nt2,tumor\nt3,tumor\n' "${id//\"/\"\"}" > "$d/metadata.csv"
  stopcase fix4_${nm}_in_sample_id "$d" 'does not match sample_id_regex'; done
# fix 5: the QC status is read BY NAME: a FAIL sample (stub hook) is left out of the merge inputs and listed, and the run still completes
F5=$(newin failqc 'a1 a2 a3 t1 t2 t3'); run "$WORK/l12e" "$WORK/12e.log" --input_dir "$F5" --mode apply --cohort development --stub_fail_sample a1; rc=$?
staged() { local f wd; f=$(ls "$1"/logs/trace_*.tsv | tail -1); wd=$(awk -F'\t' -v t="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $h["name"] ~ t {print $h["workdir"]; exit}' "$f"); ls "$wd/$3" 2>/dev/null; }
check "fix 5: FAIL sample is listed (qc_excluded.tsv) and the run completes" "[[ $rc == 0 ]] && grep -q '^a1' $F5/Result/00_input/qc_excluded.tsv"
check "fix 5: FAIL sample is NOT given to the counts merge (status read by name)" "[[ -n \"\$(staged $F5/Result MERGE_FEATURECOUNTS fc | grep a2)\" && -z \"\$(staged $F5/Result MERGE_FEATURECOUNTS fc | grep a1)\" ]]"
# fix 6: the run record holds the reference record (GENCODE release, md5 and sizes copied from manifest_quant.json)
check "fix 6: run_info.json has the reference record with the GENCODE release" "python3 -c \"import json,sys; r=json.load(open('$R1/run_info.json')).get('reference'); sys.exit(0 if r and r.get('gencode_release')=='v50' else 1)\""
check "fix 6: report manifest.json has the reference record (md5 copied)" "python3 -c \"import json,sys; m=json.load(open('$RL/manifest.json'))['reference']; sys.exit(0 if m['gencode_release']=='v50' and m['files']['gtf']['md5']=='abc123' else 1)\""
# fix 7: no duplicate MultiQC reports in 05, 06, 07; the four stage reports and the combined one exist
check "fix 7: MultiQC only in 01-04 and 08 (not in 05, 06, 07)" "[[ ! -e $R1/05_cohort_qc/qc/multiqc_report.html && ! -e $R1/06_dge/qc/multiqc_report.html && ! -e $R1/07_pathway/qc/multiqc_report.html && -e $R1/01_fastp/qc/multiqc_report.html && -e $R1/08_multiqc/multiqc_report.html ]]"
# fix 2: the report names the DESeq2 design (single batch -> no batch term)
check "fix 2: report states the DESeq2 design" "grep -q 'DESeq2 design:.*~ group.*one batch only' $RL/report/report.md"

# ---- 13. final fixes: 1 tissue wording, 3 patient-level folds, 4 sensitivity/specificity (Clopper-Pearson)
RSC=$INSTALL/env/prefix/r_classifier/bin/Rscript
# fix 1: tissue and the plain "only independent validation counts" line in the report; tissue in the index files and in DISCOVER_SUMMARY.md
check "fix 1: discover report states 'colon tumour vs normal' and the gene-selection line" "grep -q 'Intended use:.*case-vs-reference biomarker discovery' $RL/report/report.md && grep -q 'Gene selection was done on all development samples; development performance is optimistic; only independent validation counts' $RL/report/report.md"
check "fix 1: validation report has the tissue, the independence text, and NOT the development-only line" "grep -q 'case-vs-reference biomarker discovery' $WORK/report_cohort_validation/report/report.md && grep -q 'These samples were not used to build the locked model. The locked model (lock/diagnostic_v1) was applied without any change. Quality-control exclusions are listed below.' $WORK/report_cohort_validation/report/report.md && ! grep -q 'Gene selection was done' $WORK/report_cohort_validation/report/report.md"
check "fix 1: development (apply) report still has the development-only line and not the validation text" "grep -q 'Gene selection was done on all development samples; development performance is optimistic; only independent validation counts' $WORK/report_cohort_development/report/report.md && ! grep -q 'were not used to build the locked model' $WORK/report_cohort_development/report/report.md"
# report upgrade (R1-R6 and the added sections): text only, no numeric output
check "report: no empty Stages heading; Methods, Limitations present; no validation block without 04_validation" "[[ \$(awk '/^## 2. Stages/{f=1;next} f&&NF{print;exit}' $RL/report/report.md) == -* ]] && grep -q '^## Methods and versions' $RL/report/report.md && grep -q '^## Limitations and recommended next validation steps' $RL/report/report.md && grep -q 'in-sample and optimistic' $RL/report/report.md && grep -q 'Small cohort' $RL/report/report.md && ! grep -q 'Internal validation (pipeline' $RL/report/report.md && ! grep -q 'NOT built yet' $RL/report/report.md"
NM=$WORK/report_nomodel; mkdir -p $NM; cp -r $RL/input $RL/pipeline $RL/run_info.json $NM/
printf '# DISCOVER diagnostic model: NOT FITTED\nNO MODEL FITTED: only 0 candidate gene(s) passed.\n' > $NM/pipeline/discover/03_model/model_spec.txt
printf 'rule 3b: loo, 4 refits\nRULE 3B NOT EVALUABLE: only 4 refits ran (at least 10 needed); no candidate can be declared\n' > $NM/pipeline/discover/02_candidates/candidate_rules_method_note.txt
mkdir -p $NM/pipeline/discover/04_validation; printf 'metric\tvalue\tci_low\tci_high\tlabel\nAUC_nested\t0.9\t0.7\t1\tOPTIMISTIC\n' > $NM/pipeline/discover/04_validation/internal_validation_summary.tsv
printf 'WARNING: test line\n' > $NM/pipeline/quant/05_cohort_qc/cohort_qc_warnings.txt
(cd $NM && "$INSTALL/env/prefix/pyqc/bin/python3" -I "$HERE/scripts/make_report.py" --run . > "$WORK/10c.log" 2>&1)
check "report (no model): real reason (rule 3b not evaluable), no stray '(none)', no gene-selection banner" "grep -q 'NO MODEL WAS FITTED:.*RULE 3B NOT EVALUABLE' $NM/report/report.md && [[ \$(awk '/^\\*\\*NO MODEL WAS FITTED/{f=1;next} f&&NF{print;exit}' $NM/report/report.md) == '**Development performance'* ]] && ! grep -q 'Gene selection was done' $NM/report/report.md"
check "report: 04_validation block and cohort-QC WARNING lines shown only when present" "grep -q 'Internal validation (pipeline/discover/04_validation' $NM/report/report.md && grep -q 'AUC_nested' $NM/report/report.md && grep -q 'WARNING: test line' $NM/report/report.md"
check "fix 1: RESULTS_INDEX_RUN.md and the DISCOVER_SUMMARY.md template name the tissue" "grep -q 'Case vs reference groups as given in the metadata' $R1/RESULTS_INDEX_RUN.md && grep -q 'Case vs reference groups as given in the metadata' $HERE/nextflow/bin/package_model.py && grep -q 'Case vs reference groups as given in the metadata' $HERE/nextflow/main.nf"
# fix 3: with every sample its own patient the model fit is byte-identical to the previous script; with repeated patients the folds keep a patient together
FT=$WORK/fit_test; mkdir -p $FT/old $FT/new $FT/pair $FT/pair_default; git -C "$HERE" show 46c4860:nextflow/bin/fit_diagnostic_model.R > $FT/fit_old.R 2>/dev/null
(cd $FT && "$RSC" - > $FT/fixture.log 2>&1 <<'REOF'
set.seed(7); smp <- c(paste0("N", 1:6), paste0("T", 1:6)); g <- sprintf("ENSG%08d.1", 1:40)
cn <- matrix(rpois(40 * 12, 200), 40, 12, dimnames = list(g, smp)); cn[1:8, 7:12] <- cn[1:8, 7:12] * 3; cn[9:12, 1:6] <- cn[9:12, 1:6] * 3
write.table(data.frame(gene_id = g, cn, check.names = FALSE), "counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
mk <- function(pat, f) write.table(data.frame(sample_id = smp, label = smp, group = rep(c("normal", "tumour"), each = 6), baseline = rep(c("TRUE", "FALSE"), each = 6), study = "S", batch = rep(c("A", "B"), 6), count_source = "featurecounts", patient_id = pat), f, sep = "\t", quote = FALSE, row.names = FALSE)
mk(smp, "sheet_unique.tsv"); mk(rep(paste0("P", 1:6), 2), "sheet_paired.tsv")
write.table(data.frame(gene_id = g, all_rules_pass = c(rep(TRUE, 12), rep(FALSE, 28))), "cand.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_id = g, gene_name = paste0("G", 1:40), gene_type = "protein_coding", baseMean = 200, log2FoldChange_shrunk = 1), "de.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("lasso_seed: 42", "cutoff_seed: 43", "cv_folds: 3"), "q.yaml"); writeLines("routes: test", "note.txt")
REOF
)
fitrun() { (cd "$FT/$1" && "$RSC" "$2" --counts ../counts.tsv --sheet ../$3 --candidates ../cand.tsv --de ../de.tsv --questions ../q.yaml --method-note ../note.txt "${@:4}" > fit.log 2>&1); }
fitrun old "$FT/fit_old.R" sheet_unique.tsv; fitrun new "$HERE/nextflow/bin/fit_diagnostic_model.R" sheet_unique.tsv; fitrun pair "$HERE/nextflow/bin/fit_diagnostic_model.R" sheet_paired.tsv --patient-level-cv true; fitrun pair_default "$HERE/nextflow/bin/fit_diagnostic_model.R" sheet_paired.tsv
SAME=1; for f in model_genes.tsv model_spec.txt inner_fold_cutoffs.tsv development_scores.tsv development_performance.tsv; do cmp -s $FT/old/$f $FT/new/$f || SAME=0; done
if [[ -s $FT/fit_old.R ]]; then
check "fix 3: every sample its own patient -> model fit byte-identical to the previous script" "[[ $SAME == 1 && -s $FT/old/model_genes.tsv ]]"
else echo "SKIP  fix 3: byte-identical to the previous script (commit 46c4860 is not in this repository's history; not counted as pass or fail)"; fi
check "fix 3: repeated patients + --patient-level-cv true -> spec says samples of a patient are kept together (and the fit runs)" "grep -q 'all samples of one patient are kept together in the folds' $FT/pair/model_spec.txt"
check "fix 3: repeated patients, DEFAULT (switch off) -> sample-level folds, fit runs, no 'kept together' text" "[[ -s $FT/pair_default/model_genes.tsv ]] && ! grep -q 'kept together' $FT/pair_default/model_spec.txt"
(cd $FT && "$RSC" - > $FT/folds.out 2>&1 <<REOF
K <- 3; opt <- list("patient-level-cv" = "true"); pat <- rep(paste0("P", 1:6), 2); src <- readLines("$HERE/nextflow/bin/fit_diagnostic_model.R"); eval(parse(text = src[grep("^strat_foldid <- function", src):grep("^foldid <- function", src)])); REPEATED_PATIENTS <- TRUE
y <- rep(c(0, 1), each = 6); strata <- paste(y, rep(c("A", "B"), 6)); pat <- rep(paste0("P", 1:6), 2); set.seed(42); f <- foldid(strata, pat)
cat("MAXFOLDS", max(tapply(f, pat, function(v) length(unique(v)))), "NFOLDS", length(unique(f)), "\n")
REOF
)
check "fix 3: every patient's samples are in ONE fold (all K folds used)" "grep -q 'MAXFOLDS 1 NFOLDS 3' $FT/folds.out"
# fix 4: sensitivity/specificity with Clopper-Pearson CI as NEW files; the existing evaluation columns do not change
check "fix 4: development_sens_spec.tsv (discover fit) has the Clopper-Pearson columns, 6 of 6 -> lower limit 0.5407" "head -1 $FT/new/development_sens_spec.tsv | grep -q 'sensitivity.sens_ci_low.sens_ci_high' && awk -F'\t' 'NR==2 && \$3==1 && substr(\$4,1,6)==\"0.5407\"{ok=1} END{exit !ok}' $FT/new/development_sens_spec.tsv"
AP=$WORK/apply_ss; mkdir -p $AP; (cd $HERE && "$RSC" - > $AP/fixture.log 2>&1 <<REOF
gg <- read.delim("lock/diagnostic_v1/final_model_genes.tsv", stringsAsFactors = FALSE); idcol <- grep("gene_id", names(gg), value = TRUE)[1]; setcol <- names(gg)[sapply(gg, function(v) any(v %in% c("up", "down")))][1]
up <- gg[[idcol]][gg[[setcol]] == "up"]; dn <- gg[[idcol]][gg[[setcol]] == "down"]; smp <- c(paste0("n", 1:3), paste0("t", 1:3)); g <- unique(c(up, dn, sprintf("ENSG9%07d.1", 1:50)))
cn <- matrix(1000L, length(g), 6, dimnames = list(g, smp)); cn[up, 4:6] <- 40000L; cn[dn, 1:3] <- 40000L
write.table(data.frame(gene_id = g, cn, check.names = FALSE), "$AP/counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample_id = smp, group = rep(c("normal", "tumour"), each = 3), baseline = rep(c("TRUE", "FALSE"), each = 3), study = "S"), "$AP/sheet.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
REOF
)
(cd $HERE && "$RSC" scripts/apply_locked_model.R --counts $AP/counts.tsv --sample-sheet $AP/sheet.tsv --cohort development --evaluate --out-dir $AP/scores > $AP/apply.log 2>&1); rc=$?
check "fix 4: apply evaluation_sens_spec.tsv written, 3 of 3 -> sensitivity 1, lower limit 0.2924" "[[ $rc == 0 ]] && awk -F'\t' 'NR==2 && \$1==\"all\" && \$3==1 && substr(\$4,1,6)==\"0.2924\"{ok=1} END{exit !ok}' $AP/scores/evaluation_sens_spec.tsv"
check "fix 4: evaluation.tsv keeps its existing columns unchanged" "[[ \"\$(head -1 $AP/scores/evaluation.tsv)\" == \$'scope\terrors_positive\tn_positive\terrors_baseline\tn_baseline\twithin_pass_rule\tauc_not_pass_test\tlocked_pass_rule_result\tlabel' ]]"

echo; echo "$PASS passed, $FAIL failed"
exit $(( FAIL > 0 ))
