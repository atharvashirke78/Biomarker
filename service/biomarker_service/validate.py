"""Input validation = the pipeline's own scanner (scripts/input_scan.py) with the arguments the one-workflow entry passes in discover mode.
Nothing is re-implemented here, so the service and the pipeline can never disagree about what a valid input is."""
import csv
import json
import os
import re
import shutil
import subprocess
import time

from . import config


def validate_input(job_root, input_dir, mode='discover'):
    out = os.path.join(job_root, '_validate_%d' % int(time.time() * 1000))
    cmd = [config.SCANNER_PYTHON, '-I', os.path.join(config.CODE_DIR, 'scripts', 'input_scan.py'), '--input-dir', input_dir, '--mode', mode,
           '--out-dir', out, '--config', os.path.join(config.CODE_DIR, 'conf', 'input.config'), '--batch-optional', '--patient-optional',
           '--min-per-group', str(config.min_samples_per_group()), '--questions', os.path.join(config.CODE_DIR, 'config', 'questions.yaml')]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=config.CODE_DIR)
    text = (p.stdout + p.stderr).strip()
    res = dict(ok=(p.returncode == 0), exit_code=p.returncode, output=text, n_samples=None, groups=None, notes=[])
    try:
        with open(os.path.join(out, 'input_summary.json')) as f:
            s = json.load(f)
        res.update(n_samples=s.get('n_samples'), groups=s.get('groups'), notes=s.get('notes', []))
    except (OSError, ValueError):
        pass
    shutil.rmtree(out, ignore_errors=True)       # the scanner's scratch output of this check only
    if res['ok']:
        bad = unsafe_batch_values(input_dir)
        if bad:
            res['ok'] = False
            res['output'] += ('\nSTOP (service check): batch/study values may only contain letters, digits, ., _ and - (a space or comma would '
                              'make the pipeline fail after hours of computing). Please rename: ' + ', '.join(repr(b) for b in bad[:10]))
    return res


SAFE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*$')


def unsafe_batch_values(input_dir):
    """Values of batch_or_study / batch / study are put into shell commands by the pipeline (modules/sample_qc.nf) without quotes."""
    bad = []
    try:
        with open(os.path.join(input_dir, config.META_NAME), newline='', encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                for col in ('batch_or_study', 'batch', 'study'):
                    v = (row.get(col) or '').strip()
                    if v and not SAFE.match(v) and v not in bad:
                        bad.append(v)
    except OSError:
        pass
    return bad
