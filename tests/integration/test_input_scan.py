#!/usr/bin/env python3
"""Negative and positive tests of scripts/input_scan.py (folder scan, metadata.csv, group aliases, two-way check).

usage (project root): env/prefix/pyqc/bin/python3 -I tests/integration/test_input_scan.py --work <new dir>
Builds small fixture folders (tiny gzip files, not real reads) under <work>/fixtures, runs the scan per case into <work>/runs/<case>,
checks the exit code and a message fragment, and writes <work>/test_input_scan.tsv. Exit 0 only if every case passes."""
import argparse
import csv
import gzip
import os
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument('--work', required=True)
a = ap.parse_args()
if os.path.exists(a.work):
    sys.exit('ERROR: exists (never overwritten): ' + a.work)
W = os.path.realpath(a.work)
FX, RUNS = os.path.join(W, 'fixtures'), os.path.join(W, 'runs')
os.makedirs(FX)
os.makedirs(RUNS)


def fq(d, name, empty=False, plain=False):
    p = os.path.join(d, name)
    if empty:
        open(p, 'w').close()
    elif plain:
        open(p, 'w').write('@r1\nACGT\n+\nIIII\n')
    else:
        with gzip.open(p, 'wt') as f:
            f.write('@r1\nACGT\n+\nIIII\n')


def meta(d, rows, header=('sample_id', 'group', 'batch_or_study', 'patient_id')):
    with open(os.path.join(d, 'metadata.csv'), 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def folder(name, samples=('S1', 'S2', 'S3', 'S4'), tags=('_R1', '_R2'), ext='.fastq.gz'):
    d = os.path.join(FX, name)
    os.makedirs(d)
    for s in samples:
        for t in tags:
            fq(d, s + t + ext)
    return d


STD = [['S1', 'Tumor', 'B1', 'P1'], ['S2', 'Normal', 'B1', 'P1'], ['S3', 'T', 'B1', 'P2'], ['S4', 'adjacent normal', 'B1', 'P2']]
cases = []   # (name, folder, mode, expected exit, expected fragment, extra check or None)

d = folder('ok_R1R2'); meta(d, STD)
cases.append(('ok: _R1/_R2 .fastq.gz, aliases Tumor/Normal/T/adjacent normal', d, 'discover', 0, 'Tumor -> tumour (1)',
              lambda o: open(os.path.join(o, 'group_mapping.tsv')).read().count('\n') == 5))
d = folder('ok_12_fq', tags=('_1', '_2'), ext='.fq.gz'); meta(d, STD)
cases.append(('ok: _1/_2 .fq.gz', d, 'discover', 0, 'input OK: 4 samples', None))
d = folder('ok_illumina_001', tags=('_R1_001', '_R2_001')); meta(d, STD)
cases.append(('ok: Illumina _R1_001/_R2_001', d, 'discover', 0, 'input OK: 4 samples', None))
d = folder('ok_dot_R1R2', tags=('.R1', '.R2')); meta(d, STD)
cases.append(('ok: .R1/.R2', d, 'discover', 0, 'input OK: 4 samples', None))
d = folder('unpaired'); os.remove(os.path.join(d, 'S2_R2.fastq.gz')); meta(d, STD)
cases.append(('R1 without R2', d, 'discover', 2, 'S2_R1.fastq.gz: R1 without R2', None))
d = folder('mixed_styles', samples=('S1', 'S2')); fq(d, 'S3_1.fastq.gz'); fq(d, 'S3_2.fastq.gz'); meta(d, STD[:3])
cases.append(('mixed _R1/_R2 and _1/_2', d, 'discover', 2, 'mixed mate-tag styles', None))
d = folder('empty_file'); os.remove(os.path.join(d, 'S3_R1.fastq.gz')); fq(d, 'S3_R1.fastq.gz', empty=True); meta(d, STD)
cases.append(('empty file', d, 'discover', 2, 'S3_R1.fastq.gz: empty file', None))
d = folder('broken_link'); os.remove(os.path.join(d, 'S4_R2.fastq.gz')); os.symlink('/nonexistent/S4_R2.fastq.gz', os.path.join(d, 'S4_R2.fastq.gz')); meta(d, STD)
cases.append(('broken symlink', d, 'discover', 2, 'S4_R2.fastq.gz: broken symlink', None))
d = folder('not_gzip'); os.remove(os.path.join(d, 'S1_R1.fastq.gz')); fq(d, 'S1_R1.fastq.gz', plain=True); meta(d, STD)
cases.append(('not gzip', d, 'discover', 2, 'S1_R1.fastq.gz: not gzip-compressed', None))
d = folder('duplicate_id'); fq(d, 'S1_R1.fq.gz'); fq(d, 'S1_R2.fq.gz'); meta(d, STD)
cases.append(('same sample twice (.fastq.gz + .fq.gz)', d, 'discover', 2, 'given twice', None))
d = folder('no_mate_tag'); fq(d, 'S5.fastq.gz'); meta(d, STD)
cases.append(('FASTQ name without mate tag', d, 'discover', 2, 'S5.fastq.gz: no mate tag', None))
d = folder('row_without_fastq'); meta(d, STD + [['S9', 'tumour', 'B1', 'P9']])
cases.append(('metadata row without FASTQ', d, 'discover', 2, 'metadata row(s) without FASTQ: S9', None))
d = folder('fastq_without_row'); meta(d, STD[:3])
cases.append(('FASTQ without row (discover: stop)', d, 'discover', 2, 'without a metadata row: S4', None))
cases.append(('FASTQ without row (apply: unknown, listed)', d, 'apply', 0, '(empty or no metadata row) -> unknown (1)',
              lambda o: 'S4\t\tunknown' in open(os.path.join(o, 'metadata.resolved.tsv')).read()))
d = folder('no_metadata')
cases.append(('no metadata.csv in discover: draft + stop', d, 'discover', 2, 'draft written',
              lambda o: os.path.exists(os.path.join(o, 'metadata.DRAFT.csv')) and not os.path.exists(os.path.join(d, 'metadata.DRAFT.csv'))))
d = folder('unknown_label'); meta(d, STD[:3] + [['S4', 'Cancerous', 'B1', 'P2']])
cases.append(('unknown group label', d, 'discover', 2, "unrecognised group label(s): S4: 'Cancerous'", None))
d = folder('missing_column'); meta(d, [r[:2] + r[3:] for r in STD], header=('sample_id', 'group', 'patient_id'))
cases.append(('missing batch_or_study column', d, 'discover', 2, 'lacks required column(s): batch_or_study', None))
d = folder('study_override'); meta(d, [r + ['ST' + r[0]] for r in STD], header=('sample_id', 'group', 'batch_or_study', 'patient_id', 'study'))
cases.append(('separate study column overrides batch_or_study', d, 'discover', 0, 'input OK',
              lambda o: 'tumour\tFALSE\tB1\tSTS1' in open(os.path.join(o, 'metadata.resolved.tsv')).read()))
d = folder('no_patient'); meta(d, [r[:3] for r in STD], header=('sample_id', 'group', 'batch_or_study'))
cases.append(('discover without patient_id', d, 'discover', 2, 'empty patient_id', None))
d = folder('dup_rows'); meta(d, STD + [STD[0]])
cases.append(('duplicate sample_id rows', d, 'discover', 2, 'duplicate sample_id in metadata.csv: S1', None))

# ---- group aliases (case/control, disease/healthy), sample-count checks, naming styles (approved 2026-10-22 one-time improvements)
def many(n_case, n_ctrl, lab_case, lab_ctrl):
    ids = [f'C{i}' for i in range(1, n_case + 1)] + [f'N{i}' for i in range(1, n_ctrl + 1)]
    return ids, [[s, lab_case if s[0] == 'C' else lab_ctrl, 'B1', 'P' + s] for s in ids]
ids, rws = many(5, 5, 'case', 'control'); d = folder('alias_case_control', samples=ids); meta(d, rws)
cases.append(('ok: case/control -> tumour/normal', d, 'discover', 0, 'case -> tumour (5)',
              lambda o: 'control\tnormal' in open(os.path.join(o, 'metadata.resolved.tsv')).read()))
ids, rws = many(5, 5, 'Disease', 'Healthy'); d = folder('alias_disease_healthy', samples=ids); meta(d, rws)
cases.append(('ok: Disease/Healthy (case-insensitive) -> tumour/normal', d, 'discover', 0, 'Disease -> tumour (5)', None))
ids, rws = many(5, 5, 'patient', 'control'); d = folder('alias_ambiguous', samples=ids); meta(d, rws)
cases.append(("ambiguous label 'patient' is NOT mapped", d, 'discover', 2, "unrecognised group label(s): C1: 'patient'", None))
Q = ['--min-per-group', '3', '--questions', 'config/questions.yaml']      # the arguments the one-workflow entry passes (discover)
ids, rws = many(3, 3, 'tumour', 'normal'); d = folder('small_6', samples=ids); meta(d, rws)
cases.append(('6 samples: runs (exit 0), warns that rule 3b / model are NOT EVALUABLE', d, 'discover', 0, 'WARNING: 6 samples',
              lambda o: 'NOT EVALUABLE' in open(os.path.join(o, 'input_summary.json')).read(), Q))
ids, rws = many(5, 4, 'tumour', 'normal'); d = folder('small_9', samples=ids); meta(d, rws)
cases.append(('9 samples (boundary): warns, does not stop', d, 'discover', 0, 'NOT EVALUABLE', None, Q))
ids, rws = many(5, 5, 'tumour', 'normal'); d = folder('ok_10', samples=ids); meta(d, rws)
cases.append(('10 samples (boundary): no warning', d, 'discover', 0, 'input OK: 10 samples', lambda o: 'NOT EVALUABLE' not in open(os.path.join(o, 'input_summary.json')).read(), Q))
ids, rws = many(35, 35, 'tumour', 'normal'); d = folder('big_70', samples=ids); meta(d, rws)
cases.append(('70 samples with loo: warns that the discovery stage would stop, does not stop here', d, 'discover', 0, 'is above robustness_loo_max_samples = 60', None, Q))
ids, rws = many(2, 10, 'tumour', 'normal'); d = folder('group_2_of_12', samples=ids); meta(d, rws)
cases.append(('a group of 2: per-group minimum 3 still stops', d, 'discover', 2, "too few samples in group 'tumour': 2, at least 3", None, Q))
ids, rws = many(3, 3, 'tumour', 'normal'); d = folder('no_questions_6', samples=ids); meta(d, rws)
cases.append(('without --questions: no notes (previous behaviour)', d, 'discover', 0, 'input OK: 6 samples', lambda o: 'NOT EVALUABLE' not in open(os.path.join(o, 'input_summary.json')).read(), ['--min-per-group', '3']))
d = folder('ok_12_fastqgz', tags=('_1', '_2')); meta(d, STD)
cases.append(('ok: _1/_2 .fastq.gz', d, 'discover', 0, 'input OK: 4 samples', None))

# ---- U25 batch/study value rule, R8 empty-folder message, R7/H2a design warnings
J = lambda o: open(os.path.join(o, 'input_summary.json')).read()
def rowsb(batches, groups=('tumour', 'normal')):
    return [[f'S{i+1}', groups[i % 2], b, f'P{i+1}'] for i, b in enumerate(batches)]
d = folder('batch_space'); meta(d, rowsb(['Study A'] * 4))
cases.append(('batch value with a space is rejected, names value and rule', d, 'discover', 2, "'Study A' (sample S1, S2, S3, ...). Rule: letters, digits, dot, underscore and hyphen only", None))
d = folder('batch_comma'); meta(d, rowsb(['A, B'] * 4))
cases.append(('batch value with a comma is rejected', d, 'discover', 2, "'A, B'", None))
d = folder('batch_leading_dash'); meta(d, rowsb(['-B1'] * 4))
cases.append(('batch value starting with a hyphen is rejected', d, 'discover', 2, "'-B1'", None))
d = folder('batch_semicolon_quant'); meta(d, rowsb(['B1;rm'] * 4))
cases.append(('batch value with a shell character is rejected (quant mode too)', d, 'quant', 2, "regex [A-Za-z0-9][A-Za-z0-9._-]*", None))
d = folder('study_col_bad'); meta(d, [r + ['St udy'] for r in rowsb(['B1'] * 4)], header=('sample_id', 'group', 'batch_or_study', 'patient_id', 'study'))
cases.append(('separate study column with a space is rejected', d, 'discover', 2, "'St udy'", None))
d = folder('batch_good_chars'); meta(d, rowsb(['SRP.219837-a_1'] * 2 + ['B2.x'] * 2))
cases.append(('batch values with . _ - and digits are accepted', d, 'discover', 0, 'input OK: 4 samples', None))
d = os.path.join(FX, 'empty_folder'); os.makedirs(d); os.makedirs(os.path.join(d, 'sub')); fq(os.path.join(d, 'sub'), 'S1_R1.fastq.gz'); fq(d, 'S1.fastq', plain=True)
cases.append(('R8: empty folder message names the likely causes', d, 'discover', 2, 'Likely causes',
              None))
cases.append(('R8: message lists the accepted name styles and gzip and sub-folder and metadata hints', d, 'discover', 2, '_R1.fastq.gz + <sample>_R2.fastq.gz', None))
cases.append(('R8: message mentions gzip, sub-folders, metadata.csv', d, 'discover', 2, 'metadata.csv is missing', lambda o: True))
ids, rws = many(6, 6, 'tumour', 'normal'); rws = [r[:2] + ['B1' if int(r[0][1:]) % 2 else 'B2'] + r[3:] for r in rws]
d = folder('balanced_12', samples=ids); meta(d, rws)
cases.append(('R7: balanced 6+6, 2 batches: table shown, no warning', d, 'discover', 0, 'input OK: 12 samples',
              lambda o: 'batch-by-group table' in J(o) and 'WARNING' not in J(o)))
ids = ['C1', 'C2', 'C3', 'N1', 'N2', 'N3']; d = folder('one_cell', samples=ids)
meta(d, [['C1', 'tumour', 'B1', 'P1'], ['C2', 'tumour', 'B1', 'P2'], ['C3', 'tumour', 'B2', 'P3'], ['N1', 'normal', 'B1', 'P4'], ['N2', 'normal', 'B2', 'P5'], ['N3', 'normal', 'B2', 'P6']])
cases.append(('R7: one-sample cells listed (B1/normal, B2/tumour), run continues', d, 'discover', 0, 'one sample in batch-by-group cell(s) B1/normal, B2/tumour',
              lambda o: 'B1: normal=1, tumour=2' in J(o)))
ids = ['C1', 'C2', 'N1', 'N2', 'N3']; d = folder('small_groups', samples=ids)
meta(d, [['C1', 'tumour', 'B1', 'P1'], ['C2', 'tumour', 'B1', 'P2'], ['N1', 'normal', 'B1', 'P3'], ['N2', 'normal', 'B1', 'P4'], ['N3', 'normal', 'B1', 'P5']])
cases.append(('R7: single batch, smallest group 2: exploratory/underpowered warning with the benchmark text', d, 'discover', 0, 'at least 6 per group',
              lambda o: 'exploratory, underpowered' in J(o)))

rows, ok_all = [], True
for i, (name, d, mode, want, frag, extra, *xargs) in enumerate(cases, 1):
    o = os.path.join(RUNS, f'{i:02d}_{os.path.basename(d)}_{mode}')
    p = subprocess.run([sys.executable, '-I', 'scripts/input_scan.py', '--input-dir', d, '--mode', mode, '--out-dir', o] + (xargs[0] if xargs else []),
                       capture_output=True, text=True)
    msg = (p.stdout + p.stderr).strip()
    ok = p.returncode == want and frag in msg and (extra is None or extra(o))
    ok_all &= ok
    rows.append(dict(case=name, mode=mode, expected_exit=want, exit=p.returncode, expected_fragment=frag, status='PASS' if ok else 'FAIL',
                     message=msg.replace('\n', ' | ')[:400]))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  (exit {p.returncode})")
with open(os.path.join(W, 'test_input_scan.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(rows)
print(f"\n{sum(r['status'] == 'PASS' for r in rows)} of {len(rows)} PASS")
sys.exit(0 if ok_all else 1)
