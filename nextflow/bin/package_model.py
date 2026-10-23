#!/usr/bin/env python3
"""Model packaging for -entry DISCOVER (Phase 1, 2026-10-15). Writes into the current folder (the run's model/ output), never into lock/.

usage: package_model.py --model-files <files...> --inputs <name=path ...> --input-checks input_checks.json --method-note candidate_rules_method_note.txt
                        --thresholds conf/thresholds.config --questions config/questions.yaml --run-name <name>
Outputs: MODEL_MANIFEST.tsv (role, file, bytes, sha256), DISCOVER_SUMMARY.md (labels: development only, PROVISIONAL values, single-route statement).
"""
import argparse
import hashlib
import json
import os
import re


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def provisional(path):
    out = []
    for i, line in enumerate(open(path), 1):
        if 'PROVISIONAL' in line and re.search(r'[=:]\s*[0-9]', line):
            out.append(f'{os.path.basename(path)}:{i}: {line.strip()}')
    return out


ap = argparse.ArgumentParser()
ap.add_argument('--model-files', nargs='+', required=True)
ap.add_argument('--inputs', nargs='+', required=True)
ap.add_argument('--input-checks', required=True)
ap.add_argument('--method-note', required=True)
ap.add_argument('--thresholds', required=True)
ap.add_argument('--questions', required=True)
ap.add_argument('--run-name', required=True)
a = ap.parse_args()

rows = [('model', os.path.basename(p), os.path.getsize(p), sha(p)) for p in a.model_files]
for kv in a.inputs:
    k, p = kv.split('=', 1)
    rows.append((f'input:{k}', p, os.path.getsize(p), sha(p)))
with open('MODEL_MANIFEST.tsv', 'w') as o:
    o.write('role\tfile\tbytes\tsha256\n')
    for r in rows:
        o.write('\t'.join(map(str, r)) + '\n')

checks = json.load(open(a.input_checks))
note = open(a.method_note).read().strip().splitlines()
single = any('SINGLE ROUTE' in l for l in note)
spec = open('model_spec.txt').read().strip().splitlines() if os.path.exists('model_spec.txt') else []
perf = open('development_performance.tsv').read().strip().splitlines() if os.path.exists('development_performance.tsv') else []
md = [f'# DISCOVER run summary: {a.run_name}', '',
      '> **Case vs reference groups as given in the metadata. Development / technical output only. Not clinically validated.** Performance numbers below come from the training samples',
      '> (genes, rules and cutoff were chosen on them) and are OPTIMISTIC. This model is NOT the locked model `lock/diagnostic_v1`.', '']
if single:
    md += ['> **Single quantification route:** rule 2 was applied to one route only and rule 3a was SKIPPED (both need two routes).',
           '> These results are NOT equivalent to the 19-sample method, which used featureCounts and Salmon.', '']
md += ['## Input checks', f"Result: **{checks['result']}**", '']
md += [f"- [{c['status']}] {c['id']} {c['check']}: {c['detail']}" + (f" ({c['threshold']})" if c.get('threshold') else '') for c in checks['checks']]
md += ['', '## Candidate rules'] + [f'- {l}' for l in note]
# the 'candidate rules:' lines of model_spec.txt repeat the Candidate rules section above: not printed twice
md += ['', '## Model (development fit)'] + [f'    {l}' for l in spec if not l.startswith('candidate rules:')]
md += ['', '## Development performance (OPTIMISTIC, not validation)'] + [f'    {l}' for l in perf]
md += ['', '## PROVISIONAL thresholds in use (no external source)'] + [f'- {l}' for l in provisional(a.thresholds) + provisional(a.questions)]
md += ['', '## Not done in this run', '- Hold-out scoring (Phase 1 fits the model only).', '- Prognostic and predictive questions (not implemented).',
       '', '## Limitations', '- Research use only; not clinically validated.', '- Development performance is in-sample and optimistic.',
       '- Candidate genes are hypotheses and need independent validation in samples not used to build the model.',
       '', 'Every file of this model is listed with its sha256 in `MODEL_MANIFEST.tsv`.']
open('DISCOVER_SUMMARY.md', 'w').write('\n'.join(md) + '\n')
print(f'packaged {len(a.model_files)} model files and {len(a.inputs)} inputs')
