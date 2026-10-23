"""The demo: real reads (6 development samples, 100,000 read pairs each) copied into a new job, validated and queued in one click."""
import csv
import gzip
import os
import shutil

from . import config, db, validate


def available():
    return os.path.isfile(os.path.join(config.DEMO_DIR, config.META_NAME))


def info():
    if not available():
        return dict(available=False, message='demo data not found at %s' % config.DEMO_DIR)
    with open(os.path.join(config.DEMO_DIR, config.META_NAME), newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    first = sorted(n for n in os.listdir(config.DEMO_DIR) if config.FASTQ_RE.match(n))
    pairs = None
    if first:
        try:
            with gzip.open(os.path.join(config.DEMO_DIR, first[0]), 'rt') as f:
                pairs = sum(1 for _ in f) // 4
        except OSError:
            pass
    return dict(available=True, samples=[dict(sample=r['sample_id'], group=r['group'], study=r.get('batch_or_study', '')) for r in rows],
                read_pairs_per_sample=pairs, estimated_minutes=config.DEMO_MINUTES,
                note='Real reads from public colon RNA-seq studies (SRP301216, SRP219837), first 100,000 read pairs per sample. A tool demonstration: far too small for biology; '
                     'with fewer than 10 samples the biomarker model is skipped by design.')


def start():
    """Create a job from the demo files, validate it, queue it. Returns (job, validation)."""
    if not available():
        raise FileNotFoundError('demo data not found')
    jid = db.create('Demo: example sample data (6 real samples)', 'discover')
    p = db.paths(jid)
    for n in os.listdir(config.DEMO_DIR):
        if n == config.META_NAME or config.FASTQ_RE.match(n):
            shutil.copy2(os.path.join(config.DEMO_DIR, n), os.path.join(p['input'], n))
    res = validate.validate_input(p['root'], p['input'], 'discover')
    db.update(jid, validation=__import__('json').dumps(res), status='validated' if res['ok'] else 'invalid',
              message='input OK' if res['ok'] else 'demo input invalid: see the validation output')
    if res['ok']:
        db.transition(jid, ('validated',), status='queued', queued_at=db.now(), message='queued (demo)')
    return db.get(jid), res
