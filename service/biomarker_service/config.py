"""Service settings. Everything can be overridden with environment variables (used by the tests)."""
import os
import re
import sys

CODE_DIR = os.environ.get('BIOMARKER_CODE_DIR', '/workspace/Biomarker_pipeline_onewf')
INSTALL_DIR = os.environ.get('BIOMARKER_INSTALL_DIR', '/workspace/Biomarker_pipeline')
NEXTFLOW_BIN = os.environ.get('NEXTFLOW_BIN', os.path.join(INSTALL_DIR, 'env', 'bin', 'nextflow'))
SCANNER_PYTHON = os.environ.get('SERVICE_SCANNER_PYTHON', os.path.join(INSTALL_DIR, 'env', 'prefix', 'pyqc', 'bin', 'python3'))
JOBS_DIR = os.environ.get('SERVICE_JOBS_DIR', '/workspace/biomarker_jobs')
HOST = os.environ.get('SERVICE_HOST', '127.0.0.1')       # local only by default: there is no user management in this thin version
PORT = int(os.environ.get('SERVICE_PORT', '8080'))
ALLOW_REMOTE = os.environ.get('SERVICE_ALLOW_REMOTE', '0') == '1'   # there is NO login: binding to anything but localhost is refused unless this is set on purpose
API_TOKEN = os.environ.get('SERVICE_API_TOKEN', '')       # when set, every request needs it (Authorization: Bearer <token>, or the cookie set by opening /?token=<token>); required with SERVICE_ALLOW_REMOTE=1
ALLOWED_HOSTS = {'localhost', '127.0.0.1', '::1'} | {h.strip().lower() for h in os.environ.get('SERVICE_ALLOWED_HOSTS', '').split(',') if h.strip()}   # accepted Host header names (stops DNS rebinding)
RETENTION_DAYS = float(os.environ.get('SERVICE_RETENTION_DAYS', '0'))      # 0 = keep everything; N = the worker deletes finished/failed/cancelled/invalid jobs older than N days (never queued or running ones)
MAX_FILE_GB = float(os.environ.get('SERVICE_MAX_FILE_GB', '40'))
# one Nextflow run assumes it owns the whole pod: never start a second one while any other `nextflow run` is alive on this machine
CHECK_FOREIGN_RUNS = os.environ.get('SERVICE_CHECK_FOREIGN_RUNS', '1') == '1'
POLL_SECONDS = float(os.environ.get('SERVICE_POLL_SECONDS', '5'))
DEMO_DIR = os.environ.get('SERVICE_DEMO_DIR', '/workspace/biomarker_demo/input')    # real demo reads (copies of development samples), see DEMO_README.md next to it
DEMO_MINUTES = int(os.environ.get('SERVICE_DEMO_MINUTES', '10'))      # measured earlier for these 6 samples; shown to the user as an estimate
MODES = ('discover',)                                      # apply mode needs a cohort label and a locked model: not offered by the thin version

FASTQ_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,200}\.(fastq|fq)\.gz$')
META_NAME = 'metadata.csv'


def min_samples_per_group():
    """analysis_min_samples_per_group from conf/thresholds.config (the same value the pipeline's input check uses)."""
    with open(os.path.join(CODE_DIR, 'conf', 'thresholds.config')) as f:
        txt = f.read()
    m = re.search(r'^\s*analysis_min_samples_per_group\s*=\s*(\d+)', txt, re.M)
    return int(m.group(1)) if m else 3


def run_command(input_dir, output_dir, mode, resume=False):
    """The one documented pipeline command; paths are the only per-job arguments (-resume continues an interrupted job in the same folders)."""
    return [NEXTFLOW_BIN, 'run', os.path.join(CODE_DIR, 'nextflow', 'main.nf'), '--input_dir', input_dir, '--output_dir', output_dir,
            '--mode', mode, '-ansi-log', 'false'] + (['-resume'] if resume else [])
