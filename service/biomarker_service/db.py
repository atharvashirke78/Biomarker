"""Job table (sqlite). Statuses: created -> (validated | invalid) -> queued -> running -> succeeded | failed | cancelled."""
import json
import os
import sqlite3
import threading
import time
import uuid

from . import config

_lock = threading.Lock()
_job_locks = {}


def job_lock(job_id):
    """One in-process lock per job: an upload swap and the submit that queues the job cannot interleave."""
    with _lock:
        return _job_locks.setdefault(job_id, threading.Lock())


def _conn():
    os.makedirs(config.JOBS_DIR, exist_ok=True)
    c = sqlite3.connect(os.path.join(config.JOBS_DIR, 'jobs.sqlite'), timeout=30)
    c.row_factory = sqlite3.Row
    c.execute('''CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY, name TEXT, mode TEXT, status TEXT, message TEXT, created REAL, updated REAL,
                 started REAL, finished REAL, exit_code INTEGER, pid INTEGER, validation TEXT, queued_at REAL, resume INTEGER DEFAULT 0)''')
    try:
        c.execute('ALTER TABLE jobs ADD COLUMN resume INTEGER DEFAULT 0')       # databases made by an earlier version
    except sqlite3.OperationalError:
        pass
    return c


def now():
    return time.time()


def paths(job_id):
    d = os.path.join(config.JOBS_DIR, job_id)
    return dict(root=d, input=os.path.join(d, 'input'), output=os.path.join(d, 'output'), launch=os.path.join(d, 'launch'),
                console=os.path.join(d, 'console.log'))


def create(name, mode):
    job_id = uuid.uuid4().hex[:12]
    p = paths(job_id)
    os.makedirs(p['input'])
    os.makedirs(p['launch'])
    with _lock, _conn() as c:
        c.execute('INSERT INTO jobs (id, name, mode, status, message, created, updated) VALUES (?,?,?,?,?,?,?)',
                  (job_id, name, mode, 'created', 'waiting for files', now(), now()))
    return job_id


def get(job_id):
    with _lock, _conn() as c:
        r = c.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone()
    return dict(r) if r else None


def list_jobs():
    with _lock, _conn() as c:
        return [dict(r) for r in c.execute('SELECT * FROM jobs ORDER BY created DESC LIMIT 200')]


def update(job_id, **kw):
    kw['updated'] = now()
    cols = ', '.join(f'{k} = ?' for k in kw)
    with _lock, _conn() as c:
        c.execute(f'UPDATE jobs SET {cols} WHERE id = ?', (*kw.values(), job_id))


def transition(job_id, allowed_from, **kw):
    """Atomic status change: only when the current status is in allowed_from. Returns True when changed."""
    kw['updated'] = now()
    cols = ', '.join(f'{k} = ?' for k in kw)
    q = ','.join('?' for _ in allowed_from)
    with _lock, _conn() as c:
        cur = c.execute(f'UPDATE jobs SET {cols} WHERE id = ? AND status IN ({q})', (*kw.values(), job_id, *allowed_from))
        return cur.rowcount == 1


def next_queued():
    with _lock, _conn() as c:
        r = c.execute("SELECT * FROM jobs WHERE status = 'queued' ORDER BY queued_at ASC LIMIT 1").fetchone()
    return dict(r) if r else None


def running():
    with _lock, _conn() as c:
        return [dict(r) for r in c.execute("SELECT * FROM jobs WHERE status = 'running'")]


def delete(job_id):
    with _lock, _conn() as c:
        c.execute('DELETE FROM jobs WHERE id = ?', (job_id,))


def public(job):
    j = {k: job[k] for k in ('id', 'name', 'mode', 'status', 'message', 'created', 'started', 'finished', 'exit_code', 'resume')}
    j['validation'] = json.loads(job['validation']) if job.get('validation') else None
    return j


def expired(days):
    """Jobs in a final or editable state whose last update is older than `days` days (never queued or running)."""
    cutoff = now() - days * 86400
    with _lock, _conn() as c:
        return [dict(r) for r in c.execute("SELECT * FROM jobs WHERE status IN ('succeeded','failed','cancelled','invalid','created','validated') AND updated < ?", (cutoff,))]
