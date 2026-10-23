"""HTTP API + the one-page UI. Standard library only.

  POST /api/jobs                       {"name": "...", "mode": "discover"}          -> job
  PUT  /api/jobs/<id>/files/<name>     raw body = the file (FASTQ .fastq.gz/.fq.gz or metadata.csv); only while the job is not queued
  POST /api/jobs/<id>/validate         runs the pipeline's own input scan            -> {ok, output, n_samples, groups, notes}
  POST /api/jobs/<id>/submit           validates again, then queues the job
  POST /api/jobs/<id>/cancel           removes a queued job or stops a running one
  POST /api/jobs/<id>/resume           continues a failed or cancelled job with -resume (finished steps are reused)
  DELETE /api/jobs/<id>?confirm=yes    deletes a job that is not queued or running, with its uploads and results
  GET  /api/jobs, /api/jobs/<id>       status, progress, log tail
  GET  /api/jobs/<id>/results          report and key tables (when finished)
  GET  /api/jobs/<id>/summary          QC counts, candidate funnel, model, development performance (read from the output files)
  GET  /api/jobs/<id>/export.zip       tables, model files, report and run record of a finished job in one zip
  GET  /api/jobs/<id>/figures          the report's figures (PNG)
  GET  /api/jobs/<id>/live.html        live view fragment: workflow chart, per-sample matrix, pod resources, what each step does
  GET  /api/demo, POST /api/demo/start the real demo data (6 samples, 100,000 read pairs each): describe it / copy it into a new job, validate and queue it
  GET  /api/resources                  container memory, CPU and free disk
  GET  /api/jobs/<id>/files/<path>     download from the job's output folder (never work/)
"""
import hmac
import json
import mimetypes
import os
import re
import shutil
import tempfile
import threading
import uuid
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, unquote, urlparse

from . import config, db, demo, render, resources, stages, summary, validate, worker

UI = os.path.join(os.path.dirname(__file__), 'ui', 'index.html')
EDITABLE = ('created', 'validated', 'invalid')
KEY_RESULTS = ['10_report/report/report.html', '10_report/report/report.md', '10_report/manifest.json', '04_counts/sample_qc.tsv',
               '09_discover/02_candidates/candidate_table.tsv', '09_discover/02_candidates/funnel_counts.tsv', '09_discover/03_model/model_genes.tsv',
               '09_discover/03_model/model_spec.txt', '09_discover/03_model/development_performance.tsv', '04_counts/counts_featurecounts.tsv']
EXPORT_GLOBS = ['04_counts/sample_qc.tsv', '04_counts/counts_featurecounts.tsv', '04_counts/counts_salmon_tximport.tsv', 'contract/*.tsv', '06_dge/deseq2_results_all.tsv',
                '09_discover/01_de/*results_all.tsv', '09_discover/02_candidates/*', '09_discover/03_model/*.tsv', '09_discover/03_model/*.txt', '09_discover/03_model/*.md',
                '09_discover/04_validation/*', '09_locked_model/*', '00_input/*.tsv', '00_input/*.json', '10_report/report/report.md', '10_report/report/report.html',
                '10_report/manifest.json', '10_report/tool_versions.tsv', 'run_info.json', 'params.json']
RESEARCH_NOTICE = 'Research use only. Results are not clinically validated; candidate biomarkers need independent validation.'


def tail(path, n=40):
    try:
        with open(path, errors='replace') as f:
            return f.read().splitlines()[-n:]
    except OSError:
        return []


def sample_labels(job_id):
    """group of each sample (from the job's metadata.csv) and QC status (from the pipeline's sample_qc.tsv, once it exists)."""
    groups, qc = {}, {}
    try:
        import csv
        with open(os.path.join(db.paths(job_id)['input'], config.META_NAME), newline='', encoding='utf-8-sig') as f:
            for r in csv.DictReader(f):
                groups[(r.get('sample_id') or '').strip()] = (r.get('group') or '').strip()
    except OSError:
        pass
    rows = summary.read_tsv(job_id, '04_counts/sample_qc.tsv') or []
    for r in rows:
        qc[r.get('srr', '')] = r.get('status', '')
    return groups, qc


def export_files(job_id):
    """Relative paths of the tables, model files, report and run record that exist in the job's output (never work/, hidden files or anything outside output/)."""
    import glob
    out = os.path.realpath(db.paths(job_id)['output'])
    found = []
    for g in EXPORT_GLOBS:
        for f in sorted(glob.glob(os.path.join(out, *g.split('/')))):
            rel = os.path.relpath(f, out)
            if rel not in found and safe_output_path(job_id, rel):
                found.append(rel)
    return found


def safe_output_path(job_id, rel):
    out = os.path.realpath(db.paths(job_id)['output'])
    rel = unquote(rel)
    parts = [p for p in rel.split('/') if p]
    if not parts or any(p in ('..', '.') or p.startswith('.') for p in parts) or parts[0] == 'work':
        return None
    full = os.path.realpath(os.path.join(out, *parts))
    return full if full.startswith(out + os.sep) and os.path.isfile(full) else None


class Handler(BaseHTTPRequestHandler):
    server_version = 'BiomarkerService/0.1'

    def log_message(self, fmt, *a):
        pass

    # ---- access guard: Host allow-list (DNS rebinding), same-origin check on writes (cross-site requests), JSON only for bodies, optional token
    def host_name(self):
        h = (self.headers.get('Host') or '').strip().lower()
        return h[1:h.index(']')] if h.startswith('[') and ']' in h else h.rsplit(':', 1)[0] if h.count(':') == 1 else h

    def cookie_token(self):
        for part in (self.headers.get('Cookie') or '').split(';'):
            k, _, v = part.strip().partition('=')
            if k == 'svc_token':
                return v
        return ''

    def guard(self, write):
        """True when the request may go on; otherwise the answer has been sent."""
        if self.host_name() not in config.ALLOWED_HOSTS:
            self.send_json(403, dict(error='host not allowed (set SERVICE_ALLOWED_HOSTS for another name)'))
            return False
        origin = self.headers.get('Origin')
        if write and origin and urlparse(origin).netloc.lower() != (self.headers.get('Host') or '').strip().lower():
            self.send_json(403, dict(error='cross-origin request refused'))
            return False
        if write and (self.headers.get('Content-Length') or '0') not in ('0', '') and self.command in ('POST',) \
                and (self.headers.get('Content-Type') or '').split(';')[0].strip().lower() != 'application/json':
            self.send_json(415, dict(error='send JSON with Content-Type: application/json'))
            return False
        if config.API_TOKEN:
            auth = self.headers.get('Authorization') or ''
            given = auth[7:] if auth.startswith('Bearer ') else self.cookie_token()
            if not hmac.compare_digest(given.encode(), config.API_TOKEN.encode()):
                self.send_json(401, dict(error='token required: open /?token=<token> once in the browser, or send Authorization: Bearer <token>'))
                return False
        return True

    def end_headers(self):
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        super().end_headers()

    # ---- helpers
    def send_json(self, code, obj):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def body_json(self):
        n = int(self.headers.get('Content-Length') or 0)
        if n > 1_000_000:
            raise ValueError('body too large')
        raw = self.rfile.read(n) if n else b'{}'
        obj = json.loads(raw or b'{}')
        if not isinstance(obj, dict):
            raise ValueError('body must be a JSON object')
        return obj

    def route(self):
        path = urlparse(self.path).path
        return [p for p in path.split('/') if p], path

    def job_or_404(self, job_id):
        if not re.fullmatch(r'[0-9a-f]{12}', job_id or ''):
            self.send_json(404, dict(error='no such job'))
            return None
        j = db.get(job_id)
        if not j:
            self.send_json(404, dict(error='no such job'))
        return j

    def job_view(self, j):
        v = db.public(j)
        v['progress'] = worker.progress(j['id']) if j['status'] in ('running', 'succeeded', 'failed', 'cancelled') else None
        v['log_tail'] = tail(db.paths(j['id'])['console']) if j['status'] != 'created' else []
        v['input_files'] = sorted(os.listdir(db.paths(j['id'])['input']))
        v['stages'] = stages.compute((v['progress'] or {}).get('by_process', {}), (v['validation'] or {}).get('n_samples'), j['status']) if v['progress'] else None
        v['notice'] = RESEARCH_NOTICE
        return v

    # ---- verbs
    def do_HEAD(self):
        if not self.guard(False):
            return
        self.send_response(200)
        self.send_header('Content-Length', '0')
        self.end_headers()

    def do_GET(self):
        parts, path = self.route()
        q = parse_qs(urlparse(self.path).query)
        if config.API_TOKEN and path == '/' and q.get('token') and self.host_name() in config.ALLOWED_HOSTS \
                and hmac.compare_digest(q['token'][0].encode(), config.API_TOKEN.encode()):
            self.send_response(302)                          # remember the token in an HttpOnly cookie and drop it from the address bar
            self.send_header('Set-Cookie', 'svc_token=%s; HttpOnly; SameSite=Strict; Path=/' % config.API_TOKEN)
            self.send_header('Location', '/')
            self.send_header('Content-Length', '0')
            self.end_headers()
            return
        if not self.guard(False):
            return
        if path in ('/', '/index.html'):
            with open(UI, 'rb') as f:
                body = f.read()
            self.send_response(200)
            self.send_header('Content-Type', 'text/html; charset=utf-8')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        if parts == ['api', 'resources']:
            return self.send_json(200, resources.snapshot())
        if parts == ['api', 'demo']:
            return self.send_json(200, demo.info())
        if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] == 'live.html':
            j = self.job_or_404(parts[2])
            if not j:
                return
            groups, qc = sample_labels(j['id'])
            body = render.live_html(j, self.job_view(j), groups, qc).encode()
            self.send_response(200)
            self.send_header('Content-Type', 'text/html; charset=utf-8')
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(body)
            return
        if parts == ['api', 'health']:
            return self.send_json(200, dict(ok=True, running=[j['id'] for j in db.running()], retention_days=config.RETENTION_DAYS or None, max_file_gb=config.MAX_FILE_GB, notice=RESEARCH_NOTICE))
        if parts == ['api', 'jobs']:
            return self.send_json(200, [db.public(j) for j in db.list_jobs()])
        if len(parts) == 3 and parts[:2] == ['api', 'jobs']:
            j = self.job_or_404(parts[2])
            return self.send_json(200, self.job_view(j)) if j else None
        if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] == 'results':
            j = self.job_or_404(parts[2])
            if not j:
                return
            if j['status'] != 'succeeded':
                return self.send_json(409, dict(error='job has not succeeded', status=j['status']))
            found = [r for r in KEY_RESULTS if safe_output_path(j['id'], r)]
            return self.send_json(200, dict(files=found, report='10_report/report/report.html', notice=RESEARCH_NOTICE))
        if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] in ('summary', 'figures'):
            j = self.job_or_404(parts[2])
            if not j:
                return
            if parts[3] == 'figures':
                return self.send_json(200, dict(figures=summary.figures(j['id']), notice=RESEARCH_NOTICE))
            return self.send_json(200, dict(summary=summary.summary(j['id']), status=j['status'], notice=RESEARCH_NOTICE))
        if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] == 'export.zip':
            j = self.job_or_404(parts[2])
            if not j:
                return
            if j['status'] != 'succeeded':
                return self.send_json(409, dict(error='job has not succeeded', status=j['status']))
            files = export_files(j['id'])
            with tempfile.TemporaryFile() as tmp:
                with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as z:
                    for rel in files:
                        z.write(safe_output_path(j['id'], rel), arcname=rel)
                    z.writestr('README.txt', RESEARCH_NOTICE + '\nJob %s (%s). Files are copies of the pipeline outputs; 10_report/manifest.json lists their sha256.\n' % (j['id'], j['name']))
                size = tmp.tell()
                tmp.seek(0)
                self.send_response(200)
                self.send_header('Content-Type', 'application/zip')
                self.send_header('Content-Disposition', 'attachment; filename="biomarker_%s.zip"' % j['id'])
                self.send_header('Content-Length', str(size))
                self.end_headers()
                shutil.copyfileobj(tmp, self.wfile, 1 << 20)
            return
        if len(parts) >= 5 and parts[:2] == ['api', 'jobs'] and parts[3] == 'files':
            j = self.job_or_404(parts[2])
            if not j:
                return
            full = safe_output_path(j['id'], '/'.join(parts[4:]))
            if not full:
                return self.send_json(404, dict(error='no such file'))
            size = os.path.getsize(full)
            self.send_response(200)
            self.send_header('Content-Type', mimetypes.guess_type(full)[0] or 'application/octet-stream')
            self.send_header('Content-Length', str(size))
            self.send_header('Content-Security-Policy', 'sandbox allow-scripts allow-popups')     # a report runs in an opaque origin: no access to this API
            self.end_headers()
            with open(full, 'rb') as f:
                shutil.copyfileobj(f, self.wfile, 1 << 20)
            return
        self.send_json(404, dict(error='not found'))

    def do_POST(self):
        if not self.guard(True):
            return
        parts, _ = self.route()
        try:
            if parts == ['api', 'demo', 'start']:
                if not demo.available():
                    return self.send_json(409, dict(error='demo data not found on this pod'))
                job, res = demo.start()
                return self.send_json(202 if res['ok'] else 422, dict(job=db.public(job), validation=res))
            if parts == ['api', 'jobs']:
                b = self.body_json()
                mode = b.get('mode', 'discover')
                name = str(b.get('name', ''))[:80] or 'unnamed'
                if mode not in config.MODES:
                    return self.send_json(400, dict(error='mode must be one of %s' % list(config.MODES)))
                return self.send_json(201, db.public(db.get(db.create(name, mode))))
            if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] in ('validate', 'submit', 'cancel', 'resume'):
                j = self.job_or_404(parts[2])
                if not j:
                    return
                act = parts[3]
                if act == 'resume':
                    if not os.path.isdir(db.paths(j['id'])['output']):
                        return self.send_json(409, dict(error='nothing to resume: the job never started'))
                    ok = db.transition(j['id'], ('failed', 'cancelled'), status='queued', queued_at=db.now(), resume=1, exit_code=None, finished=None,
                                       message='queued to resume (finished steps are reused)')
                    return self.send_json(202 if ok else 409, db.public(db.get(j['id'])) if ok else dict(error='only a failed or cancelled job can be resumed'))
                if act == 'cancel':
                    ok = worker.cancel(j['id'])
                    return self.send_json(200 if ok else 409, dict(ok=ok, status=db.get(j['id'])['status']))
                if j['status'] not in EDITABLE:
                    return self.send_json(409, dict(error='job is %s; inputs can no longer be validated or submitted' % j['status']))
                p = db.paths(j['id'])
                res = validate.validate_input(p['root'], p['input'], j['mode'])
                db.update(j['id'], validation=json.dumps(res), status='validated' if res['ok'] else 'invalid',
                          message='input OK' if res['ok'] else 'input invalid: see the validation output')
                if act == 'validate' or not res['ok']:
                    return self.send_json(200 if act == 'validate' else 422, res)
                with db.job_lock(j['id']):              # an upload that finishes now must not slip in after the job is queued
                    queued = db.transition(j['id'], ('validated',), status='queued', queued_at=db.now(), message='queued')
                if not queued:
                    return self.send_json(409, dict(error='job could not be queued (files changed during validation?)'))
                return self.send_json(202, db.public(db.get(j['id'])))
        except (ValueError, json.JSONDecodeError) as e:
            return self.send_json(400, dict(error=str(e)))
        self.send_json(404, dict(error='not found'))

    def do_DELETE(self):
        if not self.guard(True):
            return
        parts, _ = self.route()
        if not (len(parts) == 3 and parts[:2] == ['api', 'jobs']):
            return self.send_json(404, dict(error='not found'))
        j = self.job_or_404(parts[2])
        if not j:
            return
        if parse_qs(urlparse(self.path).query).get('confirm', [''])[0] != 'yes':
            return self.send_json(400, dict(error='add ?confirm=yes: this deletes the uploads and results of the job'))
        if j['status'] in ('queued', 'running'):
            return self.send_json(409, dict(error='job is %s; cancel it first' % j['status']))
        root = os.path.realpath(db.paths(j['id'])['root'])
        if os.path.dirname(root) != os.path.realpath(config.JOBS_DIR):
            return self.send_json(500, dict(error='refusing to delete outside the jobs folder'))
        shutil.rmtree(root, ignore_errors=True)
        db.delete(j['id'])
        self.send_json(200, dict(deleted=j['id']))

    def do_PUT(self):
        if not self.guard(True):
            return
        parts, _ = self.route()
        if not (len(parts) == 5 and parts[:2] == ['api', 'jobs'] and parts[3] == 'files'):
            return self.send_json(404, dict(error='not found'))
        j = self.job_or_404(parts[2])
        if not j:
            return
        name = unquote(parts[4])
        if not (name == config.META_NAME or config.FASTQ_RE.match(name)):
            return self.send_json(400, dict(error='file name not allowed: only %s and paired FASTQ named <sample>_R1.fastq.gz style (.fastq.gz or .fq.gz, letters, digits, . _ -)' % config.META_NAME))
        if j['status'] not in EDITABLE:
            return self.send_json(409, dict(error='job is %s; files can no longer be changed' % j['status']))
        n = self.headers.get('Content-Length')
        if n is None or not n.isdigit():
            return self.send_json(411, dict(error='Content-Length required'))
        n = int(n)
        limit = int(config.MAX_FILE_GB * 1e9) if name != config.META_NAME else 5_000_000
        if n > limit:
            return self.send_json(413, dict(error='file too large (limit %d bytes)' % limit))
        p = db.paths(j['id'])
        if shutil.disk_usage(p['root']).free < n * 1.1:
            return self.send_json(507, dict(error='not enough free disk space for this file'))
        tmp = os.path.join(p['input'], '.%s.%s.part' % (name, uuid.uuid4().hex[:8]))      # own temp name: two uploads of one name cannot mix their bytes
        got = 0
        with open(tmp, 'wb') as f:
            while got < n:
                chunk = self.rfile.read(min(1 << 20, n - got))
                if not chunk:
                    break
                f.write(chunk)
                got += len(chunk)
        if got != n:
            os.remove(tmp)
            return self.send_json(400, dict(error='upload incomplete'))
        with db.job_lock(j['id']):
            if db.get(j['id'])['status'] not in EDITABLE:                  # the job was queued while this file was arriving: keep the validated input
                os.remove(tmp)
                return self.send_json(409, dict(error='job was queued while this file was uploading; the file was discarded'))
            os.replace(tmp, os.path.join(p['input'], name))
            db.update(j['id'], status='created', message='files changed; validate again', validation=None)
        self.send_json(200, dict(name=name, bytes=got))


def make_server(host=None, port=None):
    h = host or config.HOST
    if h not in ('127.0.0.1', 'localhost', '::1') and not config.ALLOW_REMOTE:
        raise SystemExit('refusing to listen on %s: use localhost (SSH tunnel) or set SERVICE_ALLOW_REMOTE=1 together with SERVICE_API_TOKEN, behind a TLS reverse proxy.' % h)
    if h not in ('127.0.0.1', 'localhost', '::1') and not config.API_TOKEN:
        raise SystemExit('refusing to listen on %s without SERVICE_API_TOKEN' % h)
    srv = ThreadingHTTPServer((host or config.HOST, config.PORT if port is None else port), Handler)
    srv.worker = worker.Worker()
    srv.worker.start()
    return srv


def main():
    srv = make_server()
    print('biomarker service on http://%s:%d  jobs in %s  (research use only)' % (config.HOST, srv.server_address[1], config.JOBS_DIR), flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.worker.stop_flag.set()


if __name__ == '__main__':
    main()
