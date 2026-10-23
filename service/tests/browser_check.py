"""Real-browser check of the one-page UI (headless Firefox via Selenium). NOT part of `unittest discover` (needs selenium, firefox, geckodriver).
usage: python3 tests/browser_check.py <scratch folder>      (starts its own service on a free port; never starts the real pipeline)
Checks: page loads and the demo table renders; the About tab opens; own-data upload through the sample-table builder (chooses files, sets groups, XHR PUT with the
Origin guard on, validate); a finished job's Summary tab shows strandedness / RNA-QC / validation; the live view shows the new chart stages; the export and report links
exist; token mode (cookie set by /?token=...) loads the same page; no JavaScript error is raised while doing so."""
import gzip
import os
import subprocess
import sys
import time

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support.ui import WebDriverWait

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
scratch = os.path.abspath(sys.argv[1])
os.makedirs(scratch, exist_ok=True)
jobs = os.path.join(scratch, 'jobs')
PY = sys.executable
FF = os.environ.get('FIREFOX_BIN', os.path.join(os.path.dirname(PY), 'firefox'))
GD = os.environ.get('GECKODRIVER_BIN', os.path.join(os.path.dirname(PY), 'geckodriver'))
results = []


def check(name, ok, extra=''):
    results.append(ok)
    print('%s  %s %s' % ('PASS' if ok else 'FAIL', name, extra), flush=True)


def start_server(port, token=''):
    env = dict(os.environ, SERVICE_JOBS_DIR=jobs, SERVICE_PORT=str(port), SERVICE_CHECK_FOREIGN_RUNS='0', NEXTFLOW_BIN='/bin/false', SERVICE_POLL_SECONDS='1',
               BIOMARKER_CODE_DIR=os.path.abspath(os.path.join(HERE, '..', '..')), SERVICE_API_TOKEN=token)
    p = subprocess.Popen([PY, '-m', 'biomarker_service.server'], cwd=os.path.join(HERE, '..'), env=env, stdout=open(os.path.join(scratch, 'svc%d.log' % port), 'w'), stderr=subprocess.STDOUT)
    for _ in range(50):
        time.sleep(0.2)
        try:
            import urllib.request
            urllib.request.urlopen('http://127.0.0.1:%d/api/health' % port, timeout=2) if not token else None
            break
        except Exception:
            pass
    return p


def make_finished_job():
    os.environ.update(SERVICE_JOBS_DIR=jobs)
    from biomarker_service import config, db
    config.JOBS_DIR = jobs
    jid = db.create('browser-check finished job', 'discover')
    o = db.paths(jid)['output']
    for d in ('05_cohort_qc', '04_counts', '09_discover/04_validation', '10_report/report', '10_report/plots'):
        os.makedirs(os.path.join(o, d), exist_ok=True)
    w = lambda rel, t: open(os.path.join(o, rel), 'w').write(t)
    w('04_counts/sample_qc.tsv', 'srr\tlabel\tstudy\tgroup\tstar_unique_pct\tfc_assigned_pct\tsalmon_mapped_pct\tfastp_q30_pct\tfastp_pass_pct\tstatus\treason\nS1\tS1\tA\ttumour\t90\t80\t90\t95\t90\tPASS\t\nS2\tS2\tA\tnormal\t40\t80\t90\t95\t90\tFAIL\tSTAR unique 40.00% < FAIL limit 50%\n')
    w('05_cohort_qc/strandedness_check.tsv', 'srr\tlabel\trseqc_forward\trseqc_reverse\trseqc_failed\trseqc_call\tsalmon_format\tsalmon_call\tfeaturecounts_s\tagree\tnote\nS1\tS1\t0.9\t0.05\t0.05\tforward\tISF\tforward\t1\tyes\t\nS2\tS2\t0.02\t0.9\t0.08\treverse\tIU\tunstranded\t0\tNO\tcheck library type\n')
    w('05_cohort_qc/rna_qc_metrics.tsv', 'srr\tlabel\texonic_pct\tintronic_pct\tdup_rate_pct\ttin_median\tbias_3p_5p\trrna_pct\tmito_pct\tflags\nS1\tS1\t98\t1\t10\t60\t1\t1\t1\t\nS2\tS2\t90\t9\t50\t20\t4\t1\t1\ttin_median low (robust z -4.0)\n')
    w('09_discover/04_validation/nested_validation_summary.tsv', 'metric\tvalue\tci_low\tci_high\tlabel\nAUC_nested_pooled\t0.9\t0.7\t1\tNESTED, rules 1-2 only; NOT independent validation\n')
    w('10_report/report/report.html', '<html><body><h1>report</h1><script>window.x=1</script></body></html>')
    db.update(jid, status='succeeded', message='finished; report written')
    return jid


def driver():
    o = Options()
    o.add_argument('-headless')
    o.binary_location = FF
    d = webdriver.Firefox(options=o, service=Service(GD))
    d.set_page_load_timeout(60)
    return d


def inject_error_log(d):
    d.execute_script("window.__errs=[];window.addEventListener('error',e=>window.__errs.push(String(e.message)));window.addEventListener('unhandledrejection',e=>window.__errs.push('rejection: '+e.reason))")


def wait(d, f, t=30):
    return WebDriverWait(d, t).until(f)


def port_free(port):
    import socket
    with socket.socket() as so:
        return so.connect_ex(('127.0.0.1', port)) != 0


port = 18200
assert port_free(port) and port_free(18201), 'ports 18200/18201 are in use (a stale test server?)'
srv = start_server(port)
d = None
try:
    jid = make_finished_job()
    d = driver()
    d.get('http://127.0.0.1:%d/' % port)
    inject_error_log(d)
    check('page loads with its title', 'Biomarker' in d.title or 'biomarker' in d.page_source.lower())
    wait(d, lambda x: x.find_element(By.ID, 'demoNote').text.strip() != 'loading demo information…')
    check('demo information is shown (table or "not available" note)', len(d.find_element(By.ID, 'demoNote').text) > 10, d.find_element(By.ID, 'demoNote').text[:60])
    d.find_element(By.ID, 'nav-about').click()
    check('About tab opens', not d.find_element(By.ID, 'page-about').get_attribute('hidden'))
    d.find_element(By.ID, 'nav-run').click()
    # open the finished job from the jobs table
    wait(d, lambda x: x.find_elements(By.CSS_SELECTOR, '[data-sel]'))
    d.find_element(By.CSS_SELECTOR, '[data-sel="%s"]' % jid).click()
    wait(d, lambda x: 'Strandedness cross-check' in x.find_element(By.ID, 'tab-summary').text, 40)
    t = d.find_element(By.ID, 'tab-summary').text
    check('summary: sample QC and flagged sample', 'Sample QC' in t and 'S2' in t)
    check('summary: strandedness disagreement shown', '1 of 2 sample(s) disagree' in t)
    check('summary: RNA-quality flag shown', 'tin_median low' in t)
    check('summary: nested validation shown with its caveat', 'NOT independent validation' in t)
    wait(d, lambda x: 'Strand + RNA QC' in x.find_element(By.ID, 'live').get_attribute('innerHTML'), 20)
    check('live view lists the new chart stage', True)
    d.set_window_size(1300, 2200)
    time.sleep(1)
    d.save_screenshot(os.path.join(scratch, 'job_view.png'))
    check('export and report links visible for a succeeded job', d.find_element(By.ID, 'exportLink').is_displayed() and d.find_element(By.ID, 'reportLink').is_displayed())
    import urllib.request
    z = urllib.request.urlopen(d.find_element(By.ID, 'exportLink').get_attribute('href'), timeout=20).read()
    check('export link returns a zip', z[:2] == b'PK')
    # own data: choose files, build the sample table, set groups, upload through XHR, validate
    fq = gzip.compress(b'@r\nACGT\n+\nIIII\n')
    files = []
    for s in ('A1', 'A2', 'A3', 'B1', 'B2', 'B3'):
        for m in ('R1', 'R2'):
            p = os.path.join(scratch, '%s_%s.fastq.gz' % (s, m))
            open(p, 'wb').write(fq)
            files.append(p)
    d.find_element(By.ID, 'files').send_keys('\n'.join(files))
    wait(d, lambda x: len(x.find_elements(By.CSS_SELECTOR, '#btable select')) == 6)
    check('sample table is built from the 12 chosen files (6 samples)', True)
    for i, sel in enumerate(d.find_elements(By.CSS_SELECTOR, '#btable select')):
        sel.send_keys('tumour' if i < 3 else 'normal')
    d.find_element(By.ID, 'upload').click()
    wait(d, lambda x: 'upload finished' in x.find_element(By.ID, 'upmsg').text or 'not accepted' in x.find_element(By.ID, 'upmsg').text, 60)
    check('upload through the page succeeds (XHR PUT passes the Origin guard)', 'upload finished' in d.find_element(By.ID, 'upmsg').text, d.find_element(By.ID, 'upmsg').text)
    d.find_element(By.ID, 'validate').click()
    wait(d, lambda x: x.find_element(By.ID, 'vout').text.strip() != '', 60)
    v = d.find_element(By.ID, 'vout').text
    check('validate shows the pipeline scanner answer', 'INPUT OK' in v or 'INPUT NOT ACCEPTED' in v, v.split('\n')[0][:80])
    errs = d.execute_script('return window.__errs')
    check('no JavaScript error was raised', not errs, str(errs))
finally:
    if d:
        d.quit()
    srv.terminate()

# token mode: the cookie set by /?token=... lets the page and its API calls through
port2 = 18201
srv2 = start_server(port2, token='s3cret')
time.sleep(2)
d = None
try:
    d = driver()
    d.get('http://127.0.0.1:%d/' % port2)
    check('token mode: the page without a token is refused', 'token required' in d.page_source)
    d.get('http://127.0.0.1:%d/?token=s3cret' % port2)
    inject_error_log(d)
    wait(d, lambda x: x.find_element(By.ID, 'demoNote').text.strip() != 'loading demo information…')
    wait(d, lambda x: x.find_elements(By.CSS_SELECTOR, '[data-sel]'))
    check('token mode: with the token the page and its API calls work (cookie)', len(d.find_elements(By.CSS_SELECTOR, '[data-sel]')) >= 1)
finally:
    if d:
        d.quit()
    srv2.terminate()
print('RESULT %d of %d checks passed' % (sum(results), len(results)))
sys.exit(0 if all(results) else 1)
