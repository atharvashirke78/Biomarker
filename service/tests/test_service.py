"""Service tests with a stand-in for Nextflow (the real pipeline is never started). Run:  <pyqc python> -m unittest discover -s service/tests -v"""
import gzip
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
import unittest
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix='svc_test_', dir=os.environ.get('SVC_TEST_TMP'))
os.environ.update(SERVICE_JOBS_DIR=os.path.join(TMP, 'jobs'), NEXTFLOW_BIN=os.path.join(HERE, 'fake_nextflow.sh'), SERVICE_SCANNER_PYTHON=sys.executable,
                  SERVICE_CHECK_FOREIGN_RUNS='0', SERVICE_POLL_SECONDS='0.3', BIOMARKER_CODE_DIR=os.environ.get('BIOMARKER_CODE_DIR', os.path.abspath(os.path.join(HERE, '..', '..'))))
sys.path.insert(0, os.path.join(HERE, '..'))
from biomarker_service import config, db, demo, render, resources, server, stages, worker  # noqa: E402


def req(method, url, data=None, headers=None):
    r = urllib.request.Request(url, data=data, method=method, headers=headers or {})
    try:
        with urllib.request.urlopen(r, timeout=60) as resp:
            body = resp.read()
            return resp.status, body
    except urllib.error.HTTPError as e:
        return e.code, e.read()


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.srv = server.make_server('127.0.0.1', 0)
        cls.base = 'http://127.0.0.1:%d' % cls.srv.server_address[1]
        threading.Thread(target=cls.srv.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.srv.worker.stop_flag.set()
        cls.srv.shutdown()

    def j(self, method, path, obj=None, **kw):
        data = json.dumps(obj).encode() if obj is not None else (b'' if method == 'POST' else None)
        code, body = req(method, self.base + path, data, dict({'Content-Type': 'application/json'}, **(kw.get('headers') or {})))
        try:
            return code, json.loads(body)
        except ValueError:
            return code, body

    def new_job(self, name='t'):
        code, job = self.j('POST', '/api/jobs', {'name': name})
        self.assertEqual(code, 201)
        return job['id']

    def put(self, jid, name, content):
        return req('PUT', '%s/api/jobs/%s/files/%s' % (self.base, jid, name), content)

    def fill(self, jid, n=6, marker=None, drop=None):
        """n samples (half tumour, half normal) with tiny gz FASTQ + metadata.csv; optional marker file (e.g. FAIL) written into the input folder."""
        fq = gzip.compress(b'@r\nACGT\n+\nIIII\n')
        rows = ['sample_id,group,batch_or_study,patient_id']
        for i in range(1, n + 1):
            rows.append('S%d,%s,B1,P%d' % (i, 'tumour' if i <= n // 2 else 'normal', i))
            for m in ('R1', 'R2'):
                if drop == (i, m):
                    continue
                self.assertEqual(self.put(jid, 'S%d_%s.fastq.gz' % (i, m), fq)[0], 200)
        self.assertEqual(self.put(jid, 'metadata.csv', ('\n'.join(rows) + '\n').encode())[0], 200)
        if marker:
            open(os.path.join(db.paths(jid)['input'], marker), 'w').close()

    def wait(self, jid, want, timeout=60):
        t0 = time.time()
        while time.time() - t0 < timeout:
            st = db.get(jid)['status']
            if st in want:
                return st
            time.sleep(0.2)
        self.fail('job %s stuck in %s, wanted %s' % (jid, db.get(jid)['status'], want))


class TestFlow(Base):
    def test_health_and_ui(self):
        code, h = self.j('GET', '/api/health')
        self.assertEqual((code, h['ok']), (200, True))
        self.assertIn('Research use only', h['notice'])
        code, body = req('GET', self.base + '/')
        self.assertEqual(code, 200)
        self.assertIn(b'Research use only', body)

    def test_full_flow_command_and_results(self):
        jid = self.new_job('flow')
        self.fill(jid, 6)
        code, v = self.j('POST', '/api/jobs/%s/validate' % jid)
        self.assertTrue(v['ok'], v)
        self.assertEqual(v['n_samples'], 6)
        self.assertTrue(any('NOT EVALUABLE' in n for n in v['notes']))      # small-cohort warning, but valid
        code, job = self.j('POST', '/api/jobs/%s/submit' % jid)
        self.assertEqual(code, 202)
        self.assertEqual(self.wait(jid, ('succeeded', 'failed')), 'succeeded')
        p = db.paths(jid)
        argv = open(os.path.join(p['output'], 'argv.txt')).read()
        main_nf = os.path.join(config.CODE_DIR, 'nextflow', 'main.nf')
        self.assertEqual(argv.strip(), 'input_dir=%s output_dir=%s extra=run %s --mode discover -ansi-log false' % (p['input'], p['output'], main_nf))
        code, res = self.j('GET', '/api/jobs/%s/results' % jid)
        self.assertEqual(code, 200)
        self.assertIn('10_report/report/report.html', res['files'])
        code, body = req('GET', '%s/api/jobs/%s/files/10_report/report/report.html' % (self.base, jid))
        self.assertEqual(code, 200)
        self.assertIn(b'report', body)
        code, body = req('GET', '%s/api/jobs/%s/files/10_report/plots/p.png' % (self.base, jid))    # relative link of the report resolves
        self.assertEqual(code, 200)
        code, view = self.j('GET', '/api/jobs/%s' % jid)
        self.assertEqual(view['progress']['completed'], 1)
        self.assertEqual(view['progress']['running'], ['SALMON_QUANT'])        # submitted in the console, not finished in the trace
        self.assertEqual(view['progress']['tasks_total'], 2)
        self.assertIn('Research use only', view['notice'])
        # a finished job's inputs can no longer be changed
        self.assertEqual(self.put(jid, 'metadata.csv', b'x')[0], 409)

    def test_invalid_input_never_queued(self):
        jid = self.new_job('missing mate')
        self.fill(jid, 6, drop=(2, 'R2'))
        code, v = self.j('POST', '/api/jobs/%s/submit' % jid)
        self.assertEqual(code, 422)
        self.assertFalse(v['ok'])
        self.assertIn('R1 without R2', v['output'])
        self.assertEqual(db.get(jid)['status'], 'invalid')
        time.sleep(1.0)
        self.assertFalse(os.path.exists(db.paths(jid)['output']))           # nothing was started

    def test_two_per_group_rule_and_small_cohort_not_blocked(self):
        jid = self.new_job('too small')
        self.fill(jid, 4)                                                   # 2 + 2: below the per-group minimum of 3
        code, v = self.j('POST', '/api/jobs/%s/validate' % jid)
        self.assertFalse(v['ok'])
        self.assertIn('too few samples in group', v['output'])
        jid2 = self.new_job('six')
        self.fill(jid2, 6)                                                  # 3 + 3: accepted (no minimum total)
        self.assertTrue(self.j('POST', '/api/jobs/%s/validate' % jid2)[1]['ok'])

    def test_upload_rules(self):
        jid = self.new_job('names')
        for bad in ('..%2Fevil.fastq.gz', 'a%2Fb.fastq.gz', 'notes.txt', '.hidden.fastq.gz', 'x.fastq'):
            self.assertEqual(self.put(jid, bad, b'abc')[0], 400, bad)
        self.assertEqual(self.put(jid, 'ok_R1.fastq.gz', b'abc')[0], 200)
        self.assertEqual(sorted(os.listdir(db.paths(jid)['input'])), ['ok_R1.fastq.gz'])
        code, _ = self.j('POST', '/api/jobs', {'name': 'x', 'mode': 'apply'})
        self.assertEqual(code, 400)                                         # apply mode is not offered

    def test_downloads_are_confined(self):
        jid = self.new_job('dl')
        self.fill(jid, 6)
        self.j('POST', '/api/jobs/%s/submit' % jid)
        self.wait(jid, ('succeeded',))
        for rel in ('work/ab/secret.txt', '../input/metadata.csv', '..%2Finput%2Fmetadata.csv', '.hidden', 'logs'):
            code, _ = req('GET', '%s/api/jobs/%s/files/%s' % (self.base, jid, rel))
            self.assertEqual(code, 404, rel)
        self.assertEqual(req('GET', self.base + '/api/jobs/zzzz/files/x')[0], 404)

    def test_failed_job_reports_the_reason(self):
        jid = self.new_job('fail')
        self.fill(jid, 6, marker='FAIL')
        self.j('POST', '/api/jobs/%s/submit' % jid)
        self.assertEqual(self.wait(jid, ('failed', 'succeeded')), 'failed')
        job = db.get(jid)
        self.assertEqual(job['exit_code'], 1)
        self.assertIn('INPUT CHECK FAILED', job['message'])

    def test_exit_zero_without_report_is_a_failure(self):
        jid = self.new_job('noreport')
        self.fill(jid, 6, marker='NOREPORT')
        self.j('POST', '/api/jobs/%s/submit' % jid)
        self.assertEqual(self.wait(jid, ('failed', 'succeeded')), 'failed')
        self.assertIn('no report', db.get(jid)['message'])


class TestQueue(Base):
    def test_one_job_at_a_time_in_order(self):
        os.environ['FAKE_SLEEP'] = '2'
        try:
            ids = []
            for k in range(3):
                jid = self.new_job('q%d' % k)
                self.fill(jid, 6)
                self.assertEqual(self.j('POST', '/api/jobs/%s/submit' % jid)[0], 202)
                ids.append(jid)
            t0 = time.time()
            while time.time() - t0 < 120 and not all(db.get(i)['status'] == 'succeeded' for i in ids):
                time.sleep(0.2)
            rows = [db.get(i) for i in ids]
            self.assertTrue(all(r['status'] == 'succeeded' for r in rows))
            starts = [r['started'] for r in rows]
            self.assertEqual(starts, sorted(starts))                              # first in, first out
            for prev, nxt in zip(rows, rows[1:]):
                self.assertLessEqual(prev['finished'], nxt['started'] + 0.5)      # never two at the same time (0.5 s: rounding of the two timestamps)
        finally:
            os.environ.pop('FAKE_SLEEP')

    def test_cancel_queued_and_running(self):
        os.environ['FAKE_SLEEP'] = '30'
        try:
            a, b = self.new_job('a'), self.new_job('b')
            for i in (a, b):
                self.fill(i, 6)
                self.j('POST', '/api/jobs/%s/submit' % i)
            self.wait(a, ('running',))
            self.assertEqual(self.j('POST', '/api/jobs/%s/cancel' % b)[0], 200)      # still queued
            self.assertEqual(db.get(b)['status'], 'cancelled')
            pid = db.get(a)['pid']
            self.assertEqual(self.j('POST', '/api/jobs/%s/cancel' % a)[0], 200)      # running
            self.wait(a, ('cancelled',), 20)
            time.sleep(0.5)
            self.assertFalse(worker.alive(pid))
        finally:
            os.environ.pop('FAKE_SLEEP')

    def test_foreign_run_guard(self):
        old = config.CHECK_FOREIGN_RUNS
        config.CHECK_FOREIGN_RUNS = True
        p = subprocess.Popen(['bash', '-c', 'sleep 20; true # nextflow run /x/nextflow/main.nf'])      # `; true` keeps bash from exec-ing sleep (which would drop the comment from the command line)
        try:
            time.sleep(0.3)
            self.assertTrue(worker.foreign_nextflow_running())
        finally:
            p.kill()
            p.wait()
            config.CHECK_FOREIGN_RUNS = old


class TestUpgrades(Base):
    def test_stage_chart_statuses(self):
        by = {'VERIFY_FASTQ': dict(done=6, active=0, failed=0), 'FASTP': dict(done=6, active=0, failed=0), 'SALMON_QUANT': dict(done=2, active=2, failed=0),
              'STAR_ALIGN': dict(done=0, active=0, failed=1)}
        st = {x['key']: x for x in stages.compute(by, 6, 'running')}
        self.assertEqual((st['verify']['status'], st['fastp']['status'], st['salmon']['status'], st['star']['status'], st['report']['status']),
                         ('done', 'done', 'running', 'failed', 'pending'))
        self.assertEqual((st['salmon']['done'], st['salmon']['expected']), (2, 6))
        self.assertEqual(st['sort']['expected'], 12)                         # two per-sample processes
        self.assertTrue(all(x['status'] == 'done' for x in stages.compute({}, 6, 'succeeded')))
        self.assertIsNone(stages.compute({}, None, 'running')[0]['expected'])   # unknown sample count: no invented total

    def test_summary_and_figures_come_from_the_output_files(self):
        jid = self.new_job('summary')
        self.fill(jid, 6)
        self.j('POST', '/api/jobs/%s/submit' % jid)
        self.wait(jid, ('succeeded',))
        code, r = self.j('GET', '/api/jobs/%s/summary' % jid)
        sm = r['summary']
        self.assertEqual(sm['qc']['counts'], {'PASS': 1, 'WARN': 1})
        self.assertEqual(sm['qc']['flagged'][0]['sample'], 'S2')
        self.assertEqual(sm['candidates']['n_pass_all_rules'], 1)
        self.assertEqual(sm['model'], dict(fitted=True, cutoff='1.25', genes=[dict(symbol='GENEA', gene_id='ENSG1', direction='up', coef='0.5')], note=None))
        self.assertIn('OPTIMISTIC', sm['performance']['rows'][0]['label'])
        self.assertTrue(any('not validation' in l for l in sm['limitations']))
        self.assertFalse(any('validated' in l and 'not' not in l and 'need' not in l for l in sm['limitations']))
        code, f = self.j('GET', '/api/jobs/%s/figures' % jid)
        caps = {x['name']: x['caption'] for x in f['figures']}
        self.assertEqual(caps['deseq2_volcano.png'], 'Volcano plot')
        code, view = self.j('GET', '/api/jobs/%s' % jid)
        self.assertTrue(all(x['status'] == 'done' for x in view['stages']))
        code, res = self.j('GET', '/api/jobs/%s/results' % jid)
        self.assertIn('04_counts/sample_qc.tsv', res['files'])                 # the real location of the sample QC table

    def test_resume_after_failure_uses_resume_flag_and_keeps_the_log(self):
        jid = self.new_job('resume')
        self.fill(jid, 6, marker='FAIL')
        self.j('POST', '/api/jobs/%s/submit' % jid)
        self.assertEqual(self.wait(jid, ('failed', 'succeeded')), 'failed')
        p = db.paths(jid)
        self.assertEqual(self.j('POST', '/api/jobs/%s/resume' % jid)[0], 202)   # queued again
        self.assertEqual(self.wait(jid, ('failed', 'succeeded')), 'failed')     # FAIL marker still there: fails again, but with -resume
        self.assertTrue(open(os.path.join(p['output'], 'argv.txt')).read().strip().endswith('-ansi-log false -resume'))
        os.remove(os.path.join(p['input'], 'FAIL'))
        self.assertEqual(self.j('POST', '/api/jobs/%s/resume' % jid)[0], 202)
        self.assertEqual(self.wait(jid, ('failed', 'succeeded')), 'succeeded')
        self.assertEqual(self.j('POST', '/api/jobs/%s/resume' % jid)[0], 409)   # a finished job cannot be resumed
        jid2 = self.new_job('never started')
        self.assertEqual(self.j('POST', '/api/jobs/%s/resume' % jid2)[0], 409)

    def test_delete_needs_confirmation_and_refuses_active_jobs(self):
        jid = self.new_job('del')
        self.fill(jid, 6)
        self.assertEqual(req('DELETE', '%s/api/jobs/%s' % (self.base, jid))[0], 400)           # no ?confirm=yes
        self.assertTrue(os.path.isdir(db.paths(jid)['root']))
        self.assertEqual(req('DELETE', '%s/api/jobs/%s?confirm=yes' % (self.base, jid))[0], 200)
        self.assertFalse(os.path.exists(db.paths(jid)['root']))
        self.assertIsNone(db.get(jid))
        os.environ['FAKE_SLEEP'] = '10'
        try:
            a = self.new_job('active')
            self.fill(a, 6)
            self.j('POST', '/api/jobs/%s/submit' % a)
            self.wait(a, ('running',))
            self.assertEqual(req('DELETE', '%s/api/jobs/%s?confirm=yes' % (self.base, a))[0], 409)
            self.j('POST', '/api/jobs/%s/cancel' % a)
            self.wait(a, ('cancelled',), 20)
        finally:
            os.environ.pop('FAKE_SLEEP')

    def test_service_check_blocks_unsafe_batch_values(self):
        jid = self.new_job('space')
        self.fill(jid, 6)
        rows = ['sample_id,group,batch_or_study'] + ['S%d,%s,Study A' % (i, 'tumour' if i <= 3 else 'normal') for i in range(1, 7)]
        self.assertEqual(self.put(jid, 'metadata.csv', ('\n'.join(rows) + '\n').encode())[0], 200)
        code, v = self.j('POST', '/api/jobs/%s/submit' % jid)
        self.assertEqual(code, 422)
        self.assertIn("not allowed", v['output'])                  # the pipeline's own scanner rejects it first (the service check is a second guard)
        self.assertIn("'Study A'", v['output'])
        self.assertEqual(db.get(jid)['status'], 'invalid')

    def test_no_login_means_localhost_only(self):
        with self.assertRaises(SystemExit):
            server.make_server('0.0.0.0', 0)
        self.assertEqual(req('GET', self.base + '/api/jobs')[0], 200)             # and no token is needed locally

    def test_export_zip_retention_and_upload_race(self):
        import io, zipfile, time
        jid = self.new_job('export')
        pth = db.paths(jid)
        os.makedirs(os.path.join(pth['output'], '09_discover', '03_model'))
        os.makedirs(os.path.join(pth['output'], 'work'))
        for rel in ('09_discover/03_model/model_genes.tsv', '04_counts_x.tsv', 'work/secret.tsv', '10_report_manifest.json'):
            os.makedirs(os.path.dirname(os.path.join(pth['output'], rel)), exist_ok=True)
            open(os.path.join(pth['output'], rel), 'w').write('x\n')
        os.makedirs(os.path.join(pth['output'], '10_report'))
        open(os.path.join(pth['output'], '10_report', 'manifest.json'), 'w').write('{}')
        self.assertEqual(req('GET', '%s/api/jobs/%s/export.zip' % (self.base, jid))[0], 409)        # not finished
        db.update(jid, status='succeeded')
        code, body = req('GET', '%s/api/jobs/%s/export.zip' % (self.base, jid))
        self.assertEqual(code, 200)
        names = zipfile.ZipFile(io.BytesIO(body)).namelist()
        self.assertIn('09_discover/03_model/model_genes.tsv', names)
        self.assertIn('10_report/manifest.json', names)
        self.assertIn('README.txt', names)
        self.assertFalse([n for n in names if n.startswith('work')])
        # upload that arrives after the job was queued is discarded
        jid2 = self.new_job('race')
        db.update(jid2, status='validated')
        with db.job_lock(jid2):
            db.update(jid2, status='queued')
        self.assertEqual(self.put(jid2, 'metadata.csv', b'sample_id,group\n')[0], 409)
        # retention: only old, non-active jobs are removed
        jid3 = self.new_job('old')
        db.update(jid3, status='failed')
        c = db._conn(); c.execute('UPDATE jobs SET updated = ? WHERE id = ?', (time.time() - 40 * 86400, jid3)); c.commit(); c.close()
        config.RETENTION_DAYS = 30
        try:
            worker.Worker().purge()
        finally:
            config.RETENTION_DAYS = 0
        self.assertIsNone(db.get(jid3))
        self.assertFalse(os.path.exists(db.paths(jid3)['root']))
        self.assertIsNotNone(db.get(jid))

    def test_summary_shows_strandedness_rnaqc_and_validation_and_chart_has_new_stages(self):
        jid = self.new_job('newparts')
        o = db.paths(jid)['output']
        for d in ('05_cohort_qc', '09_discover/04_validation'):
            os.makedirs(os.path.join(o, d))
        def w(rel, text):
            open(os.path.join(o, rel), 'w').write(text)
        w('05_cohort_qc/strandedness_check.tsv', 'srr\tlabel\trseqc_forward\trseqc_reverse\trseqc_failed\trseqc_call\tsalmon_format\tsalmon_call\tfeaturecounts_s\tagree\tnote\n'
                                                'S1\ts1\t0.9\t0.05\t0.05\tforward\tISF\tforward\t1\tyes\t\nS2\ts2\t0.02\t0.9\t0.08\treverse\tIU\tunstranded\t0\tNO\tcheck\n')
        w('05_cohort_qc/rna_qc_metrics.tsv', 'srr\tlabel\texonic_pct\tintronic_pct\tdup_rate_pct\ttin_median\tbias_3p_5p\trrna_pct\tmito_pct\tflags\nS1\ts1\t98\t1\t10\t60\t1\t1\t1\t\nS2\ts2\t90\t9\t50\t20\t4\t1\t1\ttin_median low (robust z -4.0)\n')
        w('09_discover/04_validation/nested_validation_summary.tsv', 'metric\tvalue\tci_low\tci_high\tlabel\nAUC_nested_pooled\t0.9\t0.7\t1\tNESTED, rules 1-2 only; NOT independent validation\n')
        code, v = self.j('GET', '/api/jobs/%s/summary' % jid)
        sm = v['summary']
        self.assertEqual((sm['strandedness']['n'], sm['strandedness']['n_disagree']), (2, 1))
        self.assertEqual(sm['strandedness']['rows'][0]['sample'], 'S2')
        self.assertEqual([r['sample'] for r in sm['rna_qc']['flagged']], ['S2'])
        self.assertIn('NOT independent', sm['validation']['nested'][0]['label'])
        keys = [x['key'] for x in stages.compute({}, 6, 'running')]
        self.assertIn('qcx', keys)
        self.assertIn('valid', keys)
        self.assertEqual(next(x for x in stages.compute({}, 6, 'running') if x['key'] == 'qcx')['expected'], 4 + 2 * 6)
        self.assertEqual(set(render.POS), set(keys))                                  # every stage has a box in the chart
        code, live = req('GET', '%s/api/jobs/%s/live.html' % (self.base, jid))
        self.assertEqual(code, 200)
        self.assertIn(b'Strand + RNA QC', live)

    def test_access_guard(self):
        self.assertEqual(req('GET', self.base + '/api/jobs', headers={'Host': 'evil.example'})[0], 403)              # DNS rebinding
        self.assertEqual(req('POST', self.base + '/api/jobs', b'{}', {'Origin': 'http://evil.example', 'Content-Type': 'application/json'})[0], 403)   # cross-site write
        self.assertEqual(req('POST', self.base + '/api/jobs', b'{"name":"x"}', {'Content-Type': 'text/plain'})[0], 415)
        self.assertEqual(req('POST', self.base + '/api/jobs', b'[1]', {'Content-Type': 'application/json'})[0], 400)  # not an object
        self.assertEqual(req('HEAD', self.base + '/')[0], 200)

    def test_token_and_sandboxed_report(self):
        config.API_TOKEN = 's3cret'
        try:
            self.assertEqual(req('GET', self.base + '/api/jobs')[0], 401)
            self.assertEqual(req('GET', self.base + '/api/jobs', headers={'Authorization': 'Bearer wrong'})[0], 401)
            self.assertEqual(req('GET', self.base + '/api/jobs', headers={'Authorization': 'Bearer s3cret'})[0], 200)
            self.assertEqual(req('GET', self.base + '/api/jobs', headers={'Cookie': 'svc_token=s3cret'})[0], 200)
        finally:
            config.API_TOKEN = ''
        with self.assertRaises(SystemExit):
            config.ALLOW_REMOTE = True
            try:
                server.make_server('0.0.0.0', 0)                                  # remote bind needs a token
            finally:
                config.ALLOW_REMOTE = False


class TestDemoAndLiveView(Base):
    def make_demo_dir(self, n=6):
        d = os.path.join(TMP, 'demo_input')
        os.makedirs(d, exist_ok=True)
        fq = gzip.compress(b'@r\nACGT\n+\nIIII\n')
        rows = ['sample_id,group,batch_or_study']
        for i in range(1, n + 1):
            rows.append('S%d,%s,B1' % (i, 'tumour' if i <= n // 2 else 'normal'))
            for m in ('R1', 'R2'):
                with open(os.path.join(d, 'S%d_%s.fastq.gz' % (i, m)), 'wb') as f:
                    f.write(fq)
        with open(os.path.join(d, 'metadata.csv'), 'w') as f:
            f.write('\n'.join(rows) + '\n')
        return d

    def test_demo_info_and_one_click_start(self):
        old = config.DEMO_DIR
        config.DEMO_DIR = os.path.join(TMP, 'does_not_exist')
        try:
            self.assertFalse(self.j('GET', '/api/demo')[1]['available'])
            self.assertEqual(self.j('POST', '/api/demo/start')[0], 409)               # clear answer when the demo data is missing
            config.DEMO_DIR = self.make_demo_dir()
            code, info = self.j('GET', '/api/demo')
            self.assertTrue(info['available'])
            self.assertEqual(len(info['samples']), 6)
            self.assertEqual(info['read_pairs_per_sample'], 1)
            self.assertEqual({x['group'] for x in info['samples']}, {'tumour', 'normal'})
            code, r = self.j('POST', '/api/demo/start')
            self.assertEqual(code, 202, r)
            jid = r['job']['id']
            self.assertTrue(r['validation']['ok'])
            self.assertIn('Demo', r['job']['name'])
            self.assertEqual(self.wait(jid, ('succeeded', 'failed')), 'succeeded')
            self.assertEqual(sorted(n for n in os.listdir(db.paths(jid)['input']) if n.endswith('.gz')), sorted(n for n in os.listdir(config.DEMO_DIR) if n.endswith('.gz')))
            self.assertTrue(all(os.path.isfile(os.path.join(config.DEMO_DIR, n)) for n in os.listdir(config.DEMO_DIR)))   # the demo source is only copied, never changed
        finally:
            config.DEMO_DIR = old

    def test_live_fragment_shows_chart_samples_steps_and_resources(self):
        jid = self.new_job('live')
        self.fill(jid, 6)
        self.j('POST', '/api/jobs/%s/submit' % jid)
        self.wait(jid, ('succeeded',))
        code, body = req('GET', '%s/api/jobs/%s/live.html' % (self.base, jid))
        self.assertEqual(code, 200)
        h = body.decode()
        self.assertIn('<svg', h)
        self.assertIn('st-done', h)
        self.assertIn('Align (STAR)', h)
        self.assertIn('What each step does', h)
        self.assertIn('Maps every read pair to the human genome', h)                  # plain-language description of a step
        self.assertIn('class="matrix"', h)
        self.assertIn('>S1<', h)                                                       # per-sample row from the trace
        self.assertIn('Memory in use by programs', h)
        self.assertIn('100%', h)                                                       # finished job
        self.assertEqual(req('GET', '%s/api/jobs/%s/live.html' % (self.base, 'zzzz'))[0], 404)
        code, r = self.j('GET', '/api/resources')
        self.assertEqual(code, 200)
        self.assertTrue({'mem_used_gb', 'mem_limit_gb', 'cpu_pct', 'cpus', 'disk_free_gb'} <= set(r))

    def test_render_functions_without_a_server(self):
        sts = stages.compute({'FASTP': dict(done=3, active=0, failed=0), 'STAR_ALIGN': dict(done=1, active=1, failed=0)}, 6, 'running')
        svg = render.chart_svg(sts)
        self.assertEqual(svg.count('<g class="st-'), len(stages.STAGES))
        self.assertIn('class="edge flow"', svg)                                        # the edge into the running stage is animated
        self.assertIn('3/6', svg)
        self.assertEqual(render.overall_percent(sts, 'running'), 0 if not any(x['expected'] for x in sts) else render.overall_percent(sts, 'running'))
        self.assertTrue(0 < render.overall_percent(sts, 'running') < 100)
        self.assertEqual(render.overall_percent(sts, 'succeeded'), 100)
        m = render.matrix_html({'A': {'FASTP': 'done', 'STAR_ALIGN': 'running'}, 'B<x>': {'FASTP': 'failed'}}, {'A': 'tumour'}, {'A': 'PASS'})
        self.assertIn('c-done', m)
        self.assertIn('c-running', m)
        self.assertIn('c-failed', m)
        self.assertIn('B&lt;x&gt;', m)                                                 # sample names are escaped
        self.assertIn('p-running', render.steps_html(sts))
        self.assertEqual(render.matrix_html({}, {}, {}).count('<table'), 0)
        self.assertEqual(render.elapsed_text(None, None, 'queued'), 'not started')
        r = resources.snapshot()
        self.assertIn('GB', render.resources_html(r))

    def test_every_step_has_a_description(self):
        for st in stages.STAGES:
            self.assertGreater(len(st['what']), 30, st['key'])
            self.assertTrue(st['tool'])
        self.assertEqual(len({st['key'] for st in stages.STAGES}), len(stages.STAGES))


if __name__ == '__main__':
    unittest.main()
