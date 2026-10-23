#!/usr/bin/env bash
# Starts the thin service in the background. Jobs, uploads and results go to SERVICE_JOBS_DIR (default /workspace/biomarker_jobs).
#   bash service/run_service.sh            start;   then open http://127.0.0.1:8080  (local only; use an SSH tunnel or a reverse proxy with TLS for remote use)
# Settings (environment): SERVICE_PORT, SERVICE_HOST, SERVICE_JOBS_DIR, SERVICE_API_TOKEN (required with SERVICE_ALLOW_REMOTE=1), SERVICE_ALLOWED_HOSTS, SERVICE_ALLOW_REMOTE, BIOMARKER_INSTALL_DIR, SERVICE_MAX_FILE_GB
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${BIOMARKER_INSTALL_DIR:-/workspace/Biomarker_pipeline}/env/prefix/pyqc/bin/python3"
LOGDIR="${SERVICE_JOBS_DIR:-/workspace/biomarker_jobs}"; mkdir -p "$LOGDIR"
cd "$HERE" || exit 1
setsid nohup "$PY" -m biomarker_service.server > "$LOGDIR/service.log" 2>&1 &
echo "started (pid $!), log: $LOGDIR/service.log"
