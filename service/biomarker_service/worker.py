"""One worker = one pipeline run at a time (a Nextflow run assumes it owns the whole pod: STAR 33 GB, Salmon 28 GB)."""
import csv
import glob
import os
import re
import shutil
import signal
import subprocess
import threading
import time

from . import config, db

CANCELLED = set()                 # job ids whose cancel was requested while running
REPORT = os.path.join('10_report', 'report', 'report.html')


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    try:                                              # a zombie child still answers kill(0)
        with open('/proc/%d/stat' % pid) as f:
            return f.read().rsplit(')', 1)[1].split()[0] != 'Z'
    except OSError:
        return True


def foreign_nextflow_running(own_pids=()):
    """Any `nextflow run ... main.nf` on this machine that is not one of our own runs."""
    if not config.CHECK_FOREIGN_RUNS:
        return False
    me = os.getpid()
    for d in os.listdir('/proc'):
        if not d.isdigit() or int(d) in (me, *own_pids):
            continue
        try:
            with open('/proc/%s/cmdline' % d, 'rb') as f:
                cmd = f.read().replace(b'\0', b' ').decode(errors='replace')
        except OSError:
            continue
        if re.search(r'nextflow\b.* run .*main\.nf', cmd) and 'pgrep' not in cmd:
            return True
    return False


def console_message(path, code):
    try:
        with open(path, errors='replace') as f:
            lines = f.read().splitlines()
    except OSError:
        lines = []
    for key in ('INPUT CHECK FAILED', 'STARTUP CHECKS FAILED', 'LOCK CHECK FAILED', 'ERROR ~', 'ERROR:'):
        for i, l in enumerate(lines):
            if key in l:
                return ' '.join(x.strip() for x in lines[i:i + 3])[:600]
    return 'pipeline exited with code %s; see the console log' % code


def progress(job_id):
    """Finished tasks come from the newest Nextflow trace (a trace lists a task only when it ends); tasks that were submitted but are not finished
    yet come from the 'Submitted process > NAME' lines of the console log (-ansi-log false)."""
    p = db.paths(job_id)
    traces = sorted(glob.glob(os.path.join(p['output'], 'logs', 'trace_*.tsv')), key=os.path.getmtime)
    rows = []
    if traces:
        with open(traces[-1], newline='') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
    stage = lambda n: re.sub(r' \(.*', '', n).split(':')[-1]
    finished = {r['name'] for r in rows if r['status'] in ('COMPLETED', 'CACHED', 'FAILED', 'ABORTED')}
    submitted = []
    try:
        with open(p['console'], errors='replace') as f:
            submitted = re.findall(r'Submitted process > (.+?)\s*$', f.read(), re.M)
    except OSError:
        pass
    active = [n for n in dict.fromkeys(submitted) if n not in finished]
    done = sum(r['status'] in ('COMPLETED', 'CACHED') for r in rows)
    bad = sum(r['status'] in ('FAILED', 'ABORTED') for r in rows)
    by = {}
    for r in rows:
        k = by.setdefault(stage(r['name']), dict(done=0, active=0, failed=0))
        if r['status'] in ('COMPLETED', 'CACHED'):
            k['done'] += 1
        elif r['status'] == 'FAILED':
            k['failed'] += 1
    for n in active:
        by.setdefault(stage(n), dict(done=0, active=0, failed=0))['active'] += 1
    samples = {}
    sample_of = lambda n: (re.search(r'\((.*)\)\s*$', n) or [None, None])[1]
    for r in rows:
        sid = sample_of(r['name'])
        if sid:
            samples.setdefault(sid, {})[stage(r['name'])] = 'done' if r['status'] in ('COMPLETED', 'CACHED') else ('failed' if r['status'] == 'FAILED' else samples.get(sid, {}).get(stage(r['name']), 'pending'))
    for n in active:
        sid = sample_of(n)
        if sid:
            samples.setdefault(sid, {})[stage(n)] = 'running'
    return dict(tasks_total=len(rows) + len(active), completed=done, running=sorted({stage(n) for n in active}), failed=bad,
                last_stage=stage(rows[-1]['name']) if rows else None, by_process=by, by_sample=samples)


class Worker(threading.Thread):
    daemon = True

    def __init__(self):
        super().__init__(name='pipeline-worker')
        self.stop_flag = threading.Event()
        self.proc = None

    def recover(self):
        for j in db.running():                       # a job left 'running' by a service restart
            if j['pid'] and alive(j['pid']):
                continue
            db.update(j['id'], status='failed', finished=db.now(), message='the service restarted while this job was running')

    def start_job(self, job):
        p = db.paths(job['id'])
        resume = bool(job.get('resume'))
        if os.path.exists(p['output']) and not resume:
            db.update(job['id'], status='failed', finished=db.now(), message='output folder already exists; not overwritten')
            return
        cmd = config.run_command(p['input'], p['output'], job['mode'], resume=resume)
        with open(p['console'], 'ab' if resume else 'wb') as log:
            self.proc = subprocess.Popen(cmd, cwd=p['launch'], stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        db.update(job['id'], status='running', started=db.now(), pid=self.proc.pid, message='running')

    def finish_job(self, job_id, code):
        p = db.paths(job_id)
        if job_id in CANCELLED:
            db.update(job_id, status='cancelled', finished=db.now(), exit_code=code, message='cancelled by the user')
        elif code == 0 and os.path.isfile(os.path.join(p['output'], REPORT)):
            db.update(job_id, status='succeeded', finished=db.now(), exit_code=0, message='finished; report written')
        elif code == 0:
            db.update(job_id, status='failed', finished=db.now(), exit_code=0, message='the pipeline exited with code 0 but wrote no report')
        else:
            db.update(job_id, status='failed', finished=db.now(), exit_code=code, message=console_message(p['console'], code))
        CANCELLED.discard(job_id)

    def run(self):
        self.recover()
        while not self.stop_flag.is_set():
            try:
                self.tick()
            except Exception as e:                     # never let the worker thread die
                print('worker error:', repr(e), flush=True)
            self.stop_flag.wait(config.POLL_SECONDS if self.proc is None else 1.0)

    def purge(self):
        """Retention: delete job folders and rows older than SERVICE_RETENTION_DAYS (off by default). Checks the path is a direct child of the jobs folder."""
        base = os.path.realpath(config.JOBS_DIR)
        for j in db.expired(config.RETENTION_DAYS):
            root = os.path.realpath(db.paths(j['id'])['root'])
            if os.path.dirname(root) == base:
                shutil.rmtree(root, ignore_errors=True)
                db.delete(j['id'])
                print('retention: deleted job', j['id'], flush=True)

    def tick(self):
        if config.RETENTION_DAYS > 0 and time.time() - getattr(self, 'last_purge', 0) > 3600:
            self.last_purge = time.time()
            self.purge()
        if self.proc is not None:
            code = self.proc.poll()
            if code is not None:
                jid = next((j['id'] for j in db.running() if j['pid'] == self.proc.pid), None)
                self.proc = None
                if jid:
                    self.finish_job(jid, code)
            return
        if db.running():
            return                                      # an orphan from an earlier service start is still alive
        job = db.next_queued()
        if not job or foreign_nextflow_running():
            return
        self.start_job(job)


def cancel(job_id):
    j = db.get(job_id)
    if not j:
        return False
    if db.transition(job_id, ('queued',), status='cancelled', finished=db.now(), message='cancelled before it started'):
        return True
    if j['status'] == 'running' and j['pid']:
        CANCELLED.add(job_id)
        try:
            os.killpg(j['pid'], signal.SIGTERM)
        except ProcessLookupError:
            pass
        db.update(job_id, message='cancel requested')
        return True
    return False
