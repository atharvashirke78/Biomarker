#!/usr/bin/env python3
"""Write an output manifest (path, bytes, sha256, mtime UTC) for every file under a folder, plus the git state and tool versions.
usage: env/prefix/pyqc/bin/python3 -I tests/regression/write_manifest.py <folder> <manifest.tsv>   (manifest must not exist; work/ dirs skipped)
"""
import datetime
import hashlib
import os
import subprocess
import sys

root, out = sys.argv[1], sys.argv[2]
if os.path.exists(out):
    sys.exit(f'ERROR: {out} exists (never overwritten)')
rows = []
for d, dirs, files in os.walk(root):
    dirs[:] = sorted(x for x in dirs if x != 'work' and not x.startswith('.'))
    for f in sorted(files):
        p = os.path.join(d, f)
        if os.path.abspath(p) == os.path.abspath(out) or os.path.islink(p):
            continue
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        mt = datetime.datetime.fromtimestamp(os.path.getmtime(p), datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
        rows.append((os.path.relpath(p, root), os.path.getsize(p), h.hexdigest(), mt))
git = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
dirty = len(subprocess.run(['git', 'status', '--porcelain'], capture_output=True, text=True).stdout.splitlines())
with open(out, 'w') as o:
    o.write(f'# manifest of {root} | written {datetime.datetime.now(datetime.timezone.utc):%Y-%m-%d %H:%M:%S} UTC | git {git} ({dirty} uncommitted paths) | tool versions: docs/tool_versions.tsv\n')
    o.write('path\tbytes\tsha256\tmtime_utc\n')
    for r in rows:
        o.write('\t'.join(map(str, r)) + '\n')
print(f'{len(rows)} files -> {out}')
