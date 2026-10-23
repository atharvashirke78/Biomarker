#!/usr/bin/env python3
"""Input layer of the folder entry (plan results_rc1prep_19/phase1_plan_20261016/PHASE1_PLAN.md section 7.2, owner decisions 9a).

usage (project root): env/prefix/pyqc/bin/python3 -I scripts/input_scan.py --input-dir <dir> --mode quant|discover|apply --out-dir <run>/input
                         [--counts <counts.tsv>] [--config conf/input.config] [--make-metadata]

1. Folder scan. The folder itself is listed (os.scandir + lstat, not Nextflow fromFilePairs), so broken links are seen.
   Problems (all listed, then exit 2): broken symlink, empty file, not gzip, FASTQ name without a mate tag, mixed mate-tag styles
   in one folder, R1 without R2 / R2 without R1, the same sample and mate twice (e.g. .fastq.gz and .fq.gz), sample ID not matching
   sample_id_regex. With --counts no FASTQ are needed: the samples are the counts matrix columns.
2. metadata.csv (next to the FASTQ). Required: sample_id, group, batch_or_study. batch_or_study fills batch and study; a separate
   batch or study column overrides it. Optional fastq_R1/fastq_R2 (must name the scanned files), patient_id, age (numeric), any other
   column (carried as text). Missing metadata.csv in quant/discover: metadata.DRAFT.csv is written into --out-dir and the run stops.
3. Group labels: mapped with group_aliases (case-insensitive, exact whole string); every mapping used goes to group_mapping.tsv;
   an unknown label stops the run. apply: an empty group, or a sample without a metadata row, becomes 'unknown' and is listed.
4. Two-way check: a metadata row without data stops every mode; data without a row stops quant/discover, is 'unknown' in apply.
Outputs in --out-dir: scan_files.tsv, problem_files.tsv, metadata.used.csv (byte copy), metadata.resolved.tsv, group_mapping.tsv,
samplesheet_quant.csv (QUANT sheet with fastq_1/fastq_2), clinical_discover.tsv, sample_sheet_apply.tsv, input_summary.json.
Exit 0 = input OK; 2 = input problem (listed); never writes into --input-dir."""
import argparse
import csv
import json
import os
import re
import shutil
import sys

STOP = 2
BATCH_RE = r'[A-Za-z0-9][A-Za-z0-9._-]*'      # batch / study values: they reach shell commands and the comma-split samplesheet (U25)
DEFAULT_BATCH = 'batch1'      # used only with --batch-optional


def read_config(path):
    """`key = value` lines inside the params block; quotes and // comments stripped (the same file Nextflow includes)."""
    cfg = {}
    for line in open(path):
        line = re.sub(r'\s//.*$', '', line.rstrip('\n'))
        m = re.match(r"\s*(\w+)\s*=\s*(.+?)\s*$", line)
        if m:
            v = m.group(2)
            cfg[m.group(1)] = v[1:-1] if len(v) >= 2 and v[0] == v[-1] and v[0] in "'\"" else v
    return cfg


def parse_aliases(s):
    amap = {}
    for part in s.split(';'):
        canon, labels = part.split(':', 1)
        for lab in labels.split(','):
            key = lab.strip().casefold()
            if key in amap and amap[key] != canon.strip():
                sys.exit(f'ERROR: conf/input.config group_aliases: label {lab!r} maps to two groups')
            amap[key] = canon.strip()
    return amap


def write_tsv(path, rows, cols):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def scan_folder(d, exts, styles, id_regex):
    """Returns (files, problems, pairs). pairs: sample -> {'R1': path, 'R2': path}."""
    files, problems, seen = [], [], {}
    used_styles = set()
    for e in sorted(os.scandir(d), key=lambda x: x.name):
        name = e.name
        if name.startswith('.'):
            continue
        ext = next((x for x in exts if name.endswith(x)), None)
        st = os.lstat(e.path)
        is_link = os.path.islink(e.path)
        target = os.readlink(e.path) if is_link else ''
        row = dict(file=name, kind='symlink' if is_link else ('dir' if e.is_dir(follow_symlinks=False) else 'file'),
                   link_target=target, size_bytes='', sample_id='', mate='', status='OK', reason='')
        if e.is_dir(follow_symlinks=False):          # sub-folders (e.g. the default output folder Result/) are never scanned
            row.update(status='IGNORED', reason='sub-folder: not scanned')
            files.append(row)
            continue
        if ext is None:
            row.update(status='IGNORED', reason='not a FASTQ name (' + ','.join(exts) + ')')
            files.append(row)
            continue

        def bad(reason):
            row.update(status='PROBLEM', reason=reason)
            problems.append(row)
        if is_link and not os.path.exists(e.path):
            bad('broken symlink (target missing)')
        elif not os.path.isfile(e.path):
            bad('not a regular file')
        else:
            size = os.stat(e.path).st_size
            row['size_bytes'] = size
            if size == 0:
                bad('empty file')
            else:
                with open(e.path, 'rb') as fh:
                    if fh.read(2) != b'\x1f\x8b':
                        bad('not gzip-compressed')
        stem = name[:-len(ext)]
        mate = None
        for i, (t1, t2) in enumerate(styles):
            if stem.endswith(t1):
                mate, sid = 'R1', stem[:-len(t1)]
            elif stem.endswith(t2):
                mate, sid = 'R2', stem[:-len(t2)]
            if mate:
                used_styles.add(f'{t1}/{t2}')
                break
        if not mate:
            if row['status'] == 'OK':
                bad('no mate tag (' + '; '.join(f'{a}/{b}' for a, b in styles) + ')')
            files.append(row)
            continue
        row.update(sample_id=sid, mate=mate)
        if id_regex and not re.fullmatch(id_regex, sid) and row['status'] == 'OK':
            bad(f'sample ID {sid!r} does not match sample_id_regex {id_regex!r}')
        key = (sid, mate)
        if key in seen and row['status'] == 'OK':
            bad(f'sample {sid} {mate} given twice (also {seen[key]})')
        seen.setdefault(key, name)
        files.append(row)
    if len(used_styles) > 1:
        for r in files:
            if r['mate'] and r['status'] == 'OK':
                r.update(status='PROBLEM', reason='mixed mate-tag styles in one folder: ' + ', '.join(sorted(used_styles)))
                problems.append(r)
    pairs = {}
    for r in files:
        if r['mate'] and r['status'] == 'OK':
            pairs.setdefault(r['sample_id'], {})[r['mate']] = os.path.join(d, r['file'])
    for sid, m in sorted(pairs.items()):
        for have, miss in (('R1', 'R2'), ('R2', 'R1')):
            if have in m and miss not in m:
                r = next(x for x in files if x['sample_id'] == sid and x['mate'] == have)
                r.update(status='PROBLEM', reason=f'{have} without {miss}')
                problems.append(r)
    pairs = {s: m for s, m in pairs.items() if 'R1' in m and 'R2' in m}
    return files, problems, pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-dir', required=True)
    ap.add_argument('--mode', required=True, choices=['quant', 'discover', 'apply'])
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--counts', help='counts matrix (TSV, gene_id + one column per sample): counts-level run, no FASTQ needed')
    ap.add_argument('--config', default='conf/input.config')
    ap.add_argument('--make-metadata', action='store_true', help='only write metadata.DRAFT.csv from the scanned FASTQ, then stop')
    ap.add_argument('--batch-optional', action='store_true',
                    help='batch_or_study may be missing or empty: those samples go into one batch %r (noted in input_summary.json)' % DEFAULT_BATCH)
    ap.add_argument('--patient-optional', action='store_true',
                    help='discover: patient_id may be missing or empty; such samples are treated as independent patients (patient_id = sample_id, noted in input_summary.json)')
    ap.add_argument('--min-per-group', type=int, default=0,
                    help='stop (quant/discover) if a group has fewer samples than this (0 = no check); the one-workflow entry passes analysis_min_samples_per_group')
    ap.add_argument('--questions', help='config/questions.yaml: only read to print early, NON-FATAL notes on which statistical stages the cohort size allows (discover mode)')
    a = ap.parse_args()
    d = os.path.realpath(a.input_dir)
    if not os.path.isdir(d):
        print(f'ERROR: --input-dir {a.input_dir} is not a folder'); return STOP
    if os.path.exists(a.out_dir):
        print(f'ERROR: {a.out_dir} exists (never overwritten)'); return STOP
    os.makedirs(a.out_dir)
    O = lambda f: os.path.join(a.out_dir, f)
    cfg = read_config(a.config)
    exts = [x.strip() for x in cfg['fastq_extensions'].split(',')]
    styles = [tuple(t.strip() for t in s.split(',')) for s in cfg['fastq_mate_tags'].split(';')]
    amap = parse_aliases(cfg['group_aliases'])
    baseline = cfg['baseline_group']
    msgs, summary = [], dict(input_dir=d, mode=a.mode, counts=a.counts, config=os.path.realpath(a.config))

    # ---- 1. data: FASTQ pairs, or counts columns
    files, problems, pairs = scan_folder(d, exts, styles, cfg.get('sample_id_regex', ''))
    scan_cols = ['file', 'kind', 'link_target', 'size_bytes', 'sample_id', 'mate', 'status', 'reason']
    write_tsv(O('scan_files.tsv'), files, scan_cols)
    write_tsv(O('problem_files.tsv'), problems, scan_cols)
    if a.counts:
        with open(a.counts) as f:
            data_ids = f.readline().rstrip('\n').split('\t')[1:]
        summary['data'] = f'counts matrix columns ({len(data_ids)})'
        if pairs:
            msgs.append(f'NOTE: --counts given; the {len(pairs)} FASTQ pairs in the folder are not used')
    else:
        data_ids = sorted(pairs)
        summary['data'] = f'FASTQ pairs ({len(data_ids)})'
        if problems:
            print(f'STOP: {len(problems)} problem file(s) in {d}:')
            for r in problems:
                print(f"  {r['file']}: {r['reason']}")
            print(f'(full listing: {O("problem_files.tsv")}; every scanned entry: {O("scan_files.tsv")})')
            return STOP
        if not data_ids:
            n_dirs = sum(1 for r in files if r['kind'] == 'dir')
            n_fq_like = sum(1 for r in files if re.search(r'\.(fastq|fq)(\.gz)?$', r['file']))
            print(f'STOP: no FASTQ pairs found in {d}. Likely causes:\n'
                  f'  1. wrong folder: this folder has {len(files)} entries ({n_dirs} sub-folder(s), {n_fq_like} FASTQ-like file(s)); check the path\n'
                  f'  2. unsupported file names: a pair needs one of the four accepted styles '
                  + ', '.join(f'<sample>{t1}{exts[0]} + <sample>{t2}{exts[0]}' for t1, t2 in styles)
                  + f' (extensions accepted: {", ".join(exts)})\n'
                  f'  3. .fastq / .fq files that are not gzip-compressed: only gzip files ({", ".join(exts)}) are read; compress them with gzip\n'
                  f'  4. FASTQ files inside sub-folders: only files directly in the folder are scanned, sub-folders are ignored\n'
                  f'  5. {cfg["metadata_file"]} is missing: it must sit in the same folder next to the FASTQ files '
                  f'({"found" if os.path.exists(os.path.join(d, cfg["metadata_file"])) else "NOT found"} here)\n'
                  f'  (every scanned entry: {O("scan_files.tsv")})'); return STOP
    summary['n_data_samples'] = len(data_ids)

    # ---- 2. metadata
    mpath = os.path.join(d, cfg['metadata_file'])
    if a.make_metadata or not os.path.exists(mpath):
        if a.make_metadata or a.mode in ('quant', 'discover'):
            draft = O('metadata.DRAFT.csv')
            with open(draft, 'w', newline='') as f:
                w = csv.writer(f, lineterminator='\n')
                w.writerow(['sample_id', 'group', 'batch_or_study', 'patient_id', 'fastq_R1', 'fastq_R2'])
                for s in data_ids:
                    m = pairs.get(s, {})
                    w.writerow([s, '', '', '', os.path.basename(m.get('R1', '')), os.path.basename(m.get('R2', ''))])
            why = '--make-metadata' if a.make_metadata else f'{cfg["metadata_file"]} not found in {d}'
            print(f'STOP ({why}): draft written to {draft}; fill group, batch_or_study (and patient_id for discover), '
                  f'save it as {mpath}, then start again.')
            return STOP
        meta_rows, header = [], []
        msgs.append(f'NOTE: no {cfg["metadata_file"]}; apply mode scores every sample with group = unknown')
    else:
        shutil.copyfile(mpath, O('metadata.used.csv'))
        with open(mpath, newline='', encoding='utf-8-sig') as f:
            rd = csv.DictReader(f)
            header = [h.strip() for h in (rd.fieldnames or [])]
            meta_rows = [{k.strip(): (v or '').strip() for k, v in r.items() if k is not None} for r in rd]
    errs = []
    if header:
        need = ['sample_id', 'group'] + ([] if a.batch_optional or {'batch', 'study'} <= set(header) else ['batch_or_study'])
        miss = [c for c in need if c not in header]
        if miss:
            errs.append(f'metadata.csv lacks required column(s): {", ".join(miss)} (required: sample_id, group'
                        + ('' if a.batch_optional else ', batch_or_study') + ')')
    ids = [r.get('sample_id', '') for r in meta_rows]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        errs.append(f'duplicate sample_id in metadata.csv: {", ".join(dup)}')
    if '' in ids:
        errs.append('metadata.csv has a row with an empty sample_id')
    if errs:
        print('STOP: ' + '\n      '.join(errs)); return STOP

    # ---- 3. two-way check
    meta_ids, dset = set(ids), set(data_ids)
    no_data = sorted(meta_ids - dset)
    no_row = sorted(dset - meta_ids)
    if no_data:
        errs.append(f'metadata row(s) without {"a counts column" if a.counts else "FASTQ"}: {", ".join(no_data)}')
    if no_row and a.mode != 'apply':
        errs.append(f'{"counts column(s)" if a.counts else "FASTQ sample(s)"} without a metadata row: {", ".join(no_row)} '
                    f'(mode {a.mode} needs a row per sample)')
    for r in meta_rows:
        m = pairs.get(r['sample_id'], {})
        for col, mate in (('fastq_R1', 'R1'), ('fastq_R2', 'R2')):
            if r.get(col) and not a.counts and r[col] != os.path.basename(m.get(mate, '')):
                errs.append(f"sample {r['sample_id']}: {col} = {r[col]!r} but the scanned {mate} file is {os.path.basename(m.get(mate, '')) or 'missing'!r}")
    if errs:
        print('STOP: ' + '\n      '.join(errs)); return STOP

    # ---- 4. resolve rows: group aliases, batch/study, types
    mapping, unknown, resolved, defaulted, no_patient, bad_batch = {}, [], [], [], [], {}
    extra_cols = [c for c in header if c not in ('sample_id', 'group', 'batch_or_study', 'batch', 'study', 'patient_id', 'fastq_R1', 'fastq_R2')]
    for r in meta_rows + [{'sample_id': s} for s in no_row]:
        sid, lab = r['sample_id'], r.get('group', '')
        if not lab and a.mode == 'apply':
            canon = 'unknown'
            mapping.setdefault(('(empty or no metadata row)', 'unknown'), []).append(sid)
        elif lab.casefold() in amap:
            canon = amap[lab.casefold()]
            mapping.setdefault((lab, canon), []).append(sid)
        else:
            unknown.append(f'{sid}: {lab!r}' if lab else f'{sid}: empty group')
            continue
        bos = r.get('batch_or_study', '')
        batch, study = r.get('batch') or bos, r.get('study') or bos
        if not batch or not study:
            if a.batch_optional:
                defaulted.append(sid)
                batch, study = batch or DEFAULT_BATCH, study or DEFAULT_BATCH
            else:
                if a.mode != 'apply':
                    errs.append(f'sample {sid}: empty batch/study (batch_or_study)')
                batch, study = batch or 'unknown', study or 'unknown'
        for col, val in (('batch', batch), ('study', study)):
            if val not in ('unknown', DEFAULT_BATCH) and not re.fullmatch(BATCH_RE, val):
                if sid not in bad_batch.setdefault(val, []):
                    bad_batch[val].append(sid)
        if r.get('age') and not re.fullmatch(r'\d+(\.\d+)?', r['age']):
            errs.append(f"sample {sid}: age {r['age']!r} is not a number")
        pid = r.get('patient_id', '')
        if a.mode == 'discover' and not pid and a.patient_optional:
            pid = sid
            no_patient.append(sid)
        elif a.mode == 'discover' and not pid:
            errs.append(f'sample {sid}: empty patient_id (discover needs it for the patient-level hold-out and the input checks)')
        m = pairs.get(sid, {})
        resolved.append(dict(sample_id=sid, group_original=lab, group=canon, baseline='TRUE' if canon == baseline else 'FALSE',
                             batch=batch, study=study, patient_id=pid, in_metadata='yes' if sid in meta_ids else 'no',
                             fastq_1=m.get('R1', ''), fastq_2=m.get('R2', ''), **{c: r.get(c, '') for c in extra_cols}))
    if unknown:
        known = sorted({f'{v} <- {k}' for k, v in amap.items()})
        errs.append('unrecognised group label(s): ' + '; '.join(unknown) + f' | allowed (conf/input.config group_aliases, case-insensitive): {", ".join(known)}')
    if bad_batch:
        errs.append('batch_or_study / batch / study value(s) not allowed: '
                    + '; '.join(f'{v!r} (sample {", ".join(ss[:3])}{", ..." if len(ss) > 3 else ""})' for v, ss in sorted(bad_batch.items()))
                    + f'. Rule: letters, digits, dot, underscore and hyphen only, starting with a letter or digit (regex {BATCH_RE}); '
                    'no spaces, commas or other characters, because the value is used in file names and shell commands. Example: rename \'Study A\' to \'Study_A\'.')
    groups = sorted({r['group'] for r in resolved if r['group'] != 'unknown'})
    if a.mode == 'discover' and (len(groups) != 2 or baseline not in groups):
        errs.append(f'discover needs exactly two groups including the baseline {baseline!r}; found {groups}')
    if a.mode == 'quant' and baseline not in groups:
        errs.append(f'quant needs the baseline group {baseline!r} among the samples; found {groups}')
    if a.min_per_group and a.mode != 'apply':
        for g in groups:
            n_g = sum(r['group'] == g for r in resolved)
            if n_g < a.min_per_group:
                errs.append(f'too few samples in group {g!r}: {n_g}, at least {a.min_per_group} per group are needed '
                            f'(analysis_min_samples_per_group in conf/thresholds.config, PROVISIONAL)')
    if a.mode == 'discover':       # design estimable: with more than one batch, every batch must contain both groups (otherwise group and batch cannot be separated)
        batches = sorted({r['batch'] for r in resolved})
        if len(batches) > 1:
            bad = []
            for b in batches:
                have = sorted({r['group'] for r in resolved if r['batch'] == b and r['group'] != 'unknown'})
                if len(have) < len(groups):
                    bad.append(f"batch {b!r} has only {have}")
            if bad:
                errs.append('group is confounded with batch: ' + '; '.join(bad) + f'. Every batch must contain both groups ({groups}), otherwise the group effect cannot be separated from the batch effect. '
                            'Add samples so that each batch has both groups, or put the samples into one batch (batch_or_study).')
    if errs:
        print('STOP: ' + '\n      '.join(errs)); return STOP
    if no_patient:
        msgs.append(f'NOTE: no patient_id for {len(no_patient)} of {len(resolved)} sample(s): each of them is treated as an independent patient '
                    f'(patient_id = sample_id); paired tumour/normal samples from one patient are not recognised as pairs in this run')
    if defaulted:
        msgs.append(f'NOTE: no batch_or_study for {len(defaulted)} of {len(resolved)} sample(s): all of them are placed in one batch {DEFAULT_BATCH!r}; '
                    f'batch effects cannot be separated from group effects in this run')
    if a.mode == 'discover':       # R7 / H2a: design-balance notes; information only, nothing here stops the run
        batches_all = sorted({r['batch'] for r in resolved})
        cellc = {(b, g): sum(r['batch'] == b and r['group'] == g for r in resolved) for b in batches_all for g in groups}
        table = '; '.join(f'{b}: ' + ', '.join(f'{g}={cellc[(b, g)]}' for g in groups) for b in batches_all)
        msgs.append(f'NOTE: batch-by-group table (samples per cell) -> {table}')
        one_cells = [f'{b}/{g}' for (b, g), n in cellc.items() if n == 1]
        if one_cells:
            msgs.append(f'WARNING: one sample in batch-by-group cell(s) {", ".join(one_cells)}. A cell with one sample cannot show whether that sample is typical; '
                        f'batch and group effects are hard to separate, and leave-one-out refits that drop it can be impossible. Table: {table}')
        single_b = [b for b in batches_all if sum(r['batch'] == b for r in resolved) == 1]
        if single_b and len(batches_all) > 1:
            msgs.append(f'WARNING: batch(es) {", ".join(single_b)} contain only one sample in total; a batch effect cannot be estimated from one sample.')
        if groups:
            n_min = min(sum(r['group'] == g for r in resolved) for g in groups)
            if n_min < 6:
                msgs.append(f'WARNING: the smallest group has {n_min} sample(s): exploratory, underpowered; published benchmarks recommend at least 6 per group. '
                            f'Results are hypothesis-generating only. The run continues.')
    if a.questions and a.mode == 'discover':       # per-method feasibility notes: information only, nothing here stops the run
        qt = open(a.questions).read()
        qv = lambda k, dflt: (re.search(r'^\s*' + k + r':\s*([A-Za-z0-9_.]+)', qt, re.M) or [None, dflt])[1]
        meth, n_all = qv('robustness_method', 'loo'), len(resolved)
        min_refits, loo_max = int(qv('robustness_min_refits', 10)), int(qv('robustness_loo_max_samples', 60))    # min_refits mirrors candidate_rules.R MIN_REFITS
        by_g = {g: sum(r['group'] == g for r in resolved) for g in groups}
        if meth == 'loo' and n_all < min_refits:
            msgs.append(f'WARNING: {n_all} samples ({by_g}). The leave-one-out robustness rule (rule 3b) needs at least {min_refits} refits, so candidate selection and the '
                        f'diagnostic model will be reported as NOT EVALUABLE and no gene panel is produced (fewer refits still if a batch would be left with one group). '
                        f'QC, counts, differential expression and pathway analysis still run. This is not evidence that there is no biomarker.')
        if meth == 'loo' and n_all > loo_max:
            msgs.append(f'WARNING: {n_all} samples is above robustness_loo_max_samples = {loo_max} with robustness_method = loo (config/questions.yaml): the run will fail at the discovery '
                        f'input check, after quantification, and no report is written. Changing the robustness method is a scientific decision for the owner; '
                        f'use --stop_after analysis to stop cleanly before the discovery stage.')
    resolved.sort(key=lambda r: data_ids.index(r['sample_id']))       # data order (counts column order / sorted FASTQ names)

    # ---- 5. outputs
    write_tsv(O('group_mapping.tsv'), [dict(label_in_metadata=k[0], canonical_group=k[1], n_samples=len(v), samples=','.join(v),
                                            rule='conf/input.config group_aliases (case-insensitive, whole string)') for k, v in mapping.items()],
              ['label_in_metadata', 'canonical_group', 'n_samples', 'samples', 'rule'])
    rcols = ['sample_id', 'group_original', 'group', 'baseline', 'batch', 'study', 'patient_id', 'in_metadata', 'fastq_1', 'fastq_2'] + extra_cols
    write_tsv(O('metadata.resolved.tsv'), resolved, rcols)
    if not a.counts:
        with open(O('samplesheet_quant.csv'), 'w', newline='') as f:
            w = csv.writer(f, lineterminator='\n')
            w.writerow(['sample', 'run', 'study', 'group', 'batch', 'baseline', 'fastq_1', 'fastq_2'])
            for r in resolved:
                w.writerow([r['sample_id'], r['sample_id'], r['study'], r['group'], r['batch'], r['baseline'], r['fastq_1'], r['fastq_2']])
    lab = lambda r: r.get('label') or r['sample_id']        # optional metadata column 'label' (display name), else the sample ID
    write_tsv(O('clinical_discover.tsv'), [dict(r, label=lab(r)) for r in resolved],
              ['sample_id', 'patient_id', 'group', 'baseline', 'batch', 'study', 'label'])
    write_tsv(O('sample_sheet_apply.tsv'), [dict(r, label=lab(r), count_source='featurecounts') for r in resolved],
              ['sample_id', 'label', 'group', 'baseline', 'study', 'batch', 'count_source'])
    summary.update(n_samples=len(resolved), groups={g: sum(r['group'] == g for r in resolved) for g in sorted({r['group'] for r in resolved})},
                   group_mapping=[f'{k[0]} -> {k[1]} ({len(v)})' for k, v in mapping.items()],
                   unknown_samples=[r['sample_id'] for r in resolved if r['group'] == 'unknown'], notes=msgs,
                   n_problem_files=len(problems), n_scanned_entries=len(files))
    json.dump(summary, open(O('input_summary.json'), 'w'), indent=2)
    print(f"input OK: {len(resolved)} samples {summary['groups']}; mapping: {'; '.join(summary['group_mapping'])}")
    for m in msgs:
        print(m)
    return 0


if __name__ == '__main__':
    sys.exit(main())
