#!/usr/bin/env python3
"""Run-folder bookkeeping for scripts/rnaseq.sh (run_info.json; resume guard; small glue steps). Never touches lock/ or the input folder.

usage (project root): env/prefix/pyqc/bin/python3 -I scripts/run_info.py <command> --run <run folder> [...]
  init            --mode M --input-dir D [--counts C] --cmd "..."   write run_info.json (mode, commit, uncommitted files, lock check, pod)
  resume-check    --mode M --input-dir D       exit 2 unless run_info.json has the same mode and input folder and metadata.csv is unchanged
  end             --status S --exit N          append an attempt record (start/end/status) to run_info.json
  filter-clinical --counts C                   clinical_discover.tsv -> clinical_discover.qcpass.tsv (rows present in the counts; the rest logged)
  labels-complete                              exit 0 if every sample has a known group and both groups are present (apply evaluation)"""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def sh(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout.strip()


def read(path):
    try:
        return open(path).read().strip()
    except OSError:
        return 'unreadable'


ap = argparse.ArgumentParser()
ap.add_argument('command', choices=['init', 'resume-check', 'end', 'filter-clinical', 'labels-complete'])
ap.add_argument('--run', required=True)
ap.add_argument('--mode')
ap.add_argument('--input-dir')
ap.add_argument('--counts', default='')
ap.add_argument('--cmd', default='')
ap.add_argument('--status')
ap.add_argument('--exit', type=int)
a = ap.parse_args()
J = os.path.join(a.run, 'run_info.json')

if a.command == 'init':
    if os.path.exists(J):
        sys.exit(f'ERROR: {J} exists (never overwritten)')
    _, commit = sh(['git', 'rev-parse', 'HEAD'])
    _, branch = sh(['git', 'branch', '--show-current'])
    _, status = sh(['git', 'status', '--porcelain'])
    lock_rc, _ = sh(['sha256sum', '-c', '--quiet', 'lock/diagnostic_v1/sha256sums_final.txt'])
    md = os.path.join(a.input_dir, 'metadata.csv')
    info = dict(run_id=os.path.basename(a.run), mode=a.mode, input_dir=a.input_dir, counts=a.counts or None, command=a.cmd,
                started=now(), git_commit=commit, git_branch=branch, uncommitted_files=status.splitlines(),
                lock_diagnostic_v1='7/7 OK' if lock_rc == 0 else 'MISMATCH', metadata_sha256=sha(md) if os.path.exists(md) else None,
                pod=dict(nproc=len(os.sched_getaffinity(0)), cgroup_memory_max=read('/sys/fs/cgroup/memory.max'), cpu_max=read('/sys/fs/cgroup/cpu.max')),
                label='development / technical output; not clinically validated', attempts=[dict(start=now(), kind='new')])
    json.dump(info, open(J, 'w'), indent=2)
    if lock_rc != 0:
        print('ERROR: lock/diagnostic_v1 hashes do not match sha256sums_final.txt'); sys.exit(2)

elif a.command == 'resume-check':
    info = json.load(open(J))
    errs = []
    if info['mode'] != a.mode:
        errs.append(f"run {info['run_id']} was started with --mode {info['mode']}, not {a.mode}")
    if os.path.realpath(info['input_dir']) != os.path.realpath(a.input_dir):
        errs.append(f"run {info['run_id']} was started on {info['input_dir']}, not {a.input_dir}")
    md = os.path.join(a.input_dir, 'metadata.csv')
    now_sha = sha(md) if os.path.exists(md) else None
    if now_sha != info.get('metadata_sha256'):
        errs.append('metadata.csv changed since the run started (start a new run instead)')
    if not os.path.exists(os.path.join(a.run, 'input', 'input_summary.json')):
        errs.append('the run stopped at the input check; nothing to resume (start a new run)')
    if errs:
        print('ERROR: resume refused: ' + '; '.join(errs)); sys.exit(2)
    info['attempts'].append(dict(start=now(), kind='resume'))
    json.dump(info, open(J, 'w'), indent=2)

elif a.command == 'end':
    info = json.load(open(J))
    info['attempts'][-1].update(end=now(), status=a.status, exit=a.exit)
    json.dump(info, open(J, 'w'), indent=2)

elif a.command == 'filter-clinical':
    import pandas as pd
    cl = pd.read_csv(os.path.join(a.run, 'input', 'clinical_discover.tsv'), sep='\t', dtype=str, keep_default_na=False)
    cols = open(a.counts).readline().rstrip('\n').split('\t')[1:]
    keep = cl[cl.sample_id.isin(cols)]
    drop = cl[~cl.sample_id.isin(cols)]
    keep.to_csv(os.path.join(a.run, 'input', 'clinical_discover.qcpass.tsv'), sep='\t', index=False)
    drop.assign(reason='not in the QUANT contract counts (QC FAIL)').to_csv(os.path.join(a.run, 'input', 'excluded_after_quant.tsv'), sep='\t', index=False)
    print(f'clinical rows kept {len(keep)}, excluded after QUANT {len(drop)}' + (f": {', '.join(drop.sample_id)}" if len(drop) else ''))

elif a.command == 'labels-complete':
    import pandas as pd
    r = pd.read_csv(os.path.join(a.run, 'input', 'metadata.resolved.tsv'), sep='\t', dtype=str, keep_default_na=False)
    ok = (r.group != 'unknown').all() and r.group.nunique() == 2
    sys.exit(0 if ok else 1)
