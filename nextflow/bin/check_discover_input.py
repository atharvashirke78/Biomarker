#!/usr/bin/env python3
"""Input contract of -entry DISCOVER (Phase 1, 2026-10-15). Stops the run (exit 2) on bad input; writes input_checks.json either way.

usage: check_discover_input.py --question diagnostic --counts <counts.tsv> [--counts-secondary <counts.tsv>] --clinical <clinical.tsv>
         [--holdout <ids.txt>] --symbols <gene_id_to_symbol.tsv> --min-samples N --max-unmapped-pct P [--loo-max-samples N --robustness loo|subsample]
Outputs (current folder):
  input_checks.json          every check with PASS / WARN / FAIL and its detail; PROVISIONAL thresholds are named as such
  sample_sheet.discover.tsv  training samples in the deseq2.R contract format (sample_id,label,group,baseline,study,batch,count_source,patient_id)
  counts_training.tsv        the primary counts restricted to training samples (a byte copy of the input when there is no hold-out)
  training_samples.txt, holdout_samples.txt
Clinical file (TSV), diagnostic: required sample_id, patient_id, group, baseline (TRUE for exactly one group); optional batch, study, label.
"""
import argparse
import json
import shutil
import sys

import numpy as np
import pandas as pd

checks = []


def add(cid, name, status, detail, threshold=None):
    checks.append({'id': cid, 'check': name, 'status': status, 'detail': detail, **({'threshold': threshold} if threshold else {})})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--question', required=True)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--counts-secondary')
    ap.add_argument('--clinical', required=True)
    ap.add_argument('--holdout')
    ap.add_argument('--symbols', required=True)
    ap.add_argument('--min-samples', type=int, required=True)
    ap.add_argument('--max-unmapped-pct', type=float, required=True)
    ap.add_argument('--robustness', default='loo')
    ap.add_argument('--loo-max-samples', type=int, default=60)
    a = ap.parse_args()

    # C0 question
    if a.question != 'diagnostic':
        add('C0', 'question implemented', 'FAIL', f"--question {a.question}: only 'diagnostic' is implemented (prognostic/predictive: Phase 2)")
        return finish(None)
    add('C0', 'question implemented', 'PASS', 'diagnostic')

    cn = pd.read_csv(a.counts, sep='\t', index_col=0)
    cl = pd.read_csv(a.clinical, sep='\t', dtype=str, keep_default_na=False)

    # C1 required columns and empty values
    req = ['sample_id', 'patient_id', 'group', 'baseline']
    miss = [c for c in req if c not in cl.columns]
    if miss:
        add('C1', 'required clinical columns', 'FAIL', f'missing: {miss} (diagnostic needs {req})')
        return finish(None)
    empty = {c: int((cl[c].str.strip() == '').sum()) for c in req}
    bad = {k: v for k, v in empty.items() if v}
    add('C1', 'required clinical columns', 'FAIL' if bad else 'PASS', f'empty values: {bad}' if bad else f'present and non-empty: {req}')

    # C2 sample IDs: unique, one-to-one with every matrix
    dup_cl = cl.sample_id[cl.sample_id.duplicated()].tolist()
    dup_cn = cn.columns[cn.columns.duplicated()].tolist()
    mats = {'counts': list(cn.columns)}
    if a.counts_secondary:
        mats['counts_secondary'] = list(pd.read_csv(a.counts_secondary, sep='\t', index_col=0, nrows=1).columns)
    mism = {k: {'only_in_matrix': sorted(set(v) - set(cl.sample_id)), 'only_in_clinical': sorted(set(cl.sample_id) - set(v))} for k, v in mats.items()}
    mism = {k: v for k, v in mism.items() if v['only_in_matrix'] or v['only_in_clinical']}
    ok2 = not (dup_cl or dup_cn or mism)
    add('C2', 'sample IDs match one-to-one', 'PASS' if ok2 else 'FAIL',
        f'{len(cn.columns)} samples' if ok2 else f'duplicated clinical {dup_cl}, duplicated matrix {dup_cn}, mismatch {mism}')

    # C3 counts content (primary: integer, non-negative, no NA, unique gene IDs)
    vals = cn.values
    probs = []
    if cn.index.duplicated().any():
        probs.append(f'{int(cn.index.duplicated().sum())} duplicated gene IDs')
    if pd.isna(vals).any():
        probs.append('NA values')
    else:
        if (vals < 0).any():
            probs.append('negative counts')
        if not (vals % 1 == 0).all():
            probs.append('non-integer counts (raw counts required)')
    add('C3', 'primary counts are raw integer counts', 'FAIL' if probs else 'PASS', '; '.join(probs) if probs else f'{cn.shape[0]} genes x {cn.shape[1]} samples')
    if not ok2 or bad:
        return finish(None)

    # C4 two groups, exactly one baseline group, baseline constant within a group
    cl['baseline_u'] = cl.baseline.str.upper()
    groups = sorted(cl.group.unique())
    bg = sorted(cl.group[cl.baseline_u == 'TRUE'].unique())
    incons = [g for g in groups if cl.baseline_u[cl.group == g].nunique() > 1]
    badval = sorted(set(cl.baseline_u) - {'TRUE', 'FALSE'})
    ok4 = len(groups) == 2 and len(bg) == 1 and not incons and not badval
    add('C4', 'exactly two groups, one baseline', 'PASS' if ok4 else 'FAIL',
        f'groups {groups}; baseline {bg}' + (f'; baseline differs within {incons}' if incons else '') + (f'; invalid baseline values {badval}' if badval else ''))
    if not ok4:
        return finish(None)

    # batch: batch column, else study, else one batch
    if 'batch' not in cl.columns or (cl.batch.str.strip() == '').all():
        cl['batch'] = cl['study'] if 'study' in cl.columns else 'batch1'
        add('W1', 'batch column', 'WARN', "no 'batch' column: batch = study" if 'study' in cl.columns else "no 'batch' or 'study' column: one batch (design ~ group)")
    if 'study' not in cl.columns:
        cl['study'] = cl['batch']
    if 'label' not in cl.columns:
        cl['label'] = cl['sample_id']

    # C5 hold-out: IDs known, patient-level separation
    hold = []
    if a.holdout:
        hold = [l.strip() for l in open(a.holdout) if l.strip()]
        unknown = sorted(set(hold) - set(cl.sample_id))
        if a.counts_secondary:
            add('C5', 'hold-out supported for this input', 'FAIL', 'a hold-out with two quantification routes is not supported in v1 (the tximport object cannot be subset here)')
            return finish(None)
        tr_pat = set(cl.patient_id[~cl.sample_id.isin(hold)]); ho_pat = set(cl.patient_id[cl.sample_id.isin(hold)])
        both = sorted(tr_pat & ho_pat)
        ok5 = not unknown and not both and len(hold) > 0
        add('C5', 'hold-out patient-level', 'PASS' if ok5 else 'FAIL',
            f'{len(hold)} hold-out samples, {len(ho_pat)} patients' if ok5 else f'unknown IDs {unknown}; patients in both training and hold-out {both}')
        if not ok5:
            return finish(None)
    else:
        add('C5', 'hold-out patient-level', 'PASS', 'no hold-out given: all samples are training; every performance number is labelled development only, optimistic')
    tr = cl[~cl.sample_id.isin(hold)].copy()

    # C6 minimum samples per group (training)
    n_by = tr.group.value_counts().to_dict()
    ok6 = all(n_by.get(g, 0) >= a.min_samples for g in groups)
    add('C6', 'minimum samples per group (training)', 'PASS' if ok6 else 'FAIL', f'{n_by}', threshold=f'>= {a.min_samples} per group (conf/thresholds.config analysis_min_samples_per_group; PROVISIONAL, no source)')

    # C7 group confounded with batch (training): every batch must contain both groups; H2b gate: the DESeq2 model matrix (~ batch + group, or ~ group
    # with one batch) must also have full column rank with residual degrees of freedom. Only turns the error DESeq2 would raise later into an early stop.
    tab = pd.crosstab(tr.batch, tr.group)
    conf = tab.index[(tab == 0).any(axis=1)].tolist()
    tab_txt = 'batch-by-group table (rows = batch, columns = group):\n' + tab.to_string()
    mm_rank = mm_cols = None
    if not conf:
        bt = sorted(tr.batch.unique()); gb = sorted(tr.group.unique(), key=lambda g: (g not in bg, g))
        mm = [[1.0] + [float(b == x) for x in bt[1:]] + [float(g == gb[1])] for b, g in zip(tr.batch, tr.group)]
        mm_rank, mm_cols = int(np.linalg.matrix_rank(np.array(mm))), len(mm[0])
    deficient = mm_rank is not None and (mm_rank < mm_cols or len(tr) - mm_rank < 1)
    if conf:
        detail = (f'batches with only one group: {conf}. {tab_txt}\nGroup is nested in batch, so ~ batch + group cannot be fitted. To change: add samples of the missing group to that batch, '
                  'give the same batch label to batches that are one experiment, or put one label in the batch column for all samples (design ~ group)')
    elif deficient:
        detail = (f'model matrix has rank {mm_rank} of {mm_cols} columns with {len(tr)} training samples ({len(tr) - mm_rank} residual degrees of freedom). {tab_txt}\n'
                  'To change: add samples, or merge batches so there are fewer batch terms than the samples allow')
    else:
        detail = f'{tab.shape[0]} batch(es), each with both groups'
    add('C7', 'group not confounded with batch', 'FAIL' if conf or deficient else 'PASS', detail)

    # C8 unmapped gene IDs (versioned match, else unversioned)
    sym = pd.read_csv(a.symbols, sep='\t')
    ids = set(sym.iloc[:, 0]); ids_u = {i.split('.')[0] for i in ids}
    unm = [g for g in cn.index if g not in ids and str(g).split('.')[0] not in ids_u]
    pct = 100 * len(unm) / len(cn.index)
    add('C8', 'gene IDs map to symbols', 'PASS' if pct <= a.max_unmapped_pct else 'FAIL', f'{pct:.2f} % unmapped ({len(unm)} of {len(cn.index)})',
        threshold=f'<= {a.max_unmapped_pct} % (conf/thresholds.config analysis_max_unmapped_gene_pct; PROVISIONAL, no source)')

    # C9 robustness method feasible
    n_tr = len(tr)
    ok9 = not (a.robustness == 'loo' and n_tr > a.loo_max_samples)
    add('C9', 'robustness method feasible', 'PASS' if ok9 else 'FAIL', f'{a.robustness} with {n_tr} training samples',
        threshold=f'loo only up to {a.loo_max_samples} samples (config/questions.yaml robustness_loo_max_samples; PROVISIONAL)')

    # warnings
    rep = tr.patient_id.value_counts(); rep = rep[rep > 1]
    if len(rep):
        add('W2', 'patients with more than one training sample', 'WARN', f'{len(rep)} patients ({int(rep.sum())} samples): samples are not independent; the model has no patient term')
    if tr.batch.nunique() == 1:
        add('W3', 'single batch', 'WARN', 'design ~ group')
    if a.counts_secondary is None:
        add('W4', 'single quantification route', 'WARN', 'rule 2 applied to one route only; rule 3a SKIPPED (config/questions.yaml); results are NOT equivalent to the 19-sample method')

    if any(c['status'] == 'FAIL' for c in checks):
        return finish(None)
    order = list(cn.columns)
    tr = tr.set_index('sample_id').loc[[s for s in order if s not in hold]].reset_index()
    tr['count_source'] = 'featurecounts'
    tr[['sample_id', 'label', 'group', 'baseline', 'study', 'batch', 'count_source', 'patient_id']].to_csv('sample_sheet.discover.tsv', sep='\t', index=False)
    if hold:
        cn[[s for s in order if s not in hold]].to_csv('counts_training.tsv', sep='\t')
    else:
        shutil.copyfile(a.counts, 'counts_training.tsv')
    open('training_samples.txt', 'w').write('\n'.join(tr.sample_id) + '\n')
    open('holdout_samples.txt', 'w').write('\n'.join(hold) + ('\n' if hold else ''))
    return finish(True)


def finish(ok):
    fail = [c for c in checks if c['status'] == 'FAIL']
    json.dump({'result': 'PASS' if ok and not fail else 'FAIL', 'checks': checks}, open('input_checks.json', 'w'), indent=2)
    for c in checks:
        print(f"[{c['status']}] {c['id']} {c['check']}: {c['detail']}" + (f" | {c['threshold']}" if c.get('threshold') else ''))
    if fail or not ok:
        print('INPUT CHECKS FAILED: the run stops here', file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
