"""Live pod resources for the demo view: container memory, CPU and free disk (read from the cgroup; nothing is measured by the pipeline itself)."""
import os
import shutil
import threading
import time

from . import config

_lock = threading.Lock()
_last = dict(t=0.0, usage=None, res=None)


def _read(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return None


def snapshot():
    """Memory used by programs (anonymous memory), file cache, limit, CPU % since the previous call (cached for 1 s), free disk."""
    with _lock:
        now = time.time()
        if _last['res'] is not None and now - _last['t'] < 1.0:
            return _last['res']
        stat = {}
        for line in (_read('/sys/fs/cgroup/memory.stat') or '').splitlines():
            k, _, v = line.partition(' ')
            if v.isdigit():
                stat[k] = int(v)
        limit = (_read('/sys/fs/cgroup/memory.max') or '').strip()
        cpu = {}
        for line in (_read('/sys/fs/cgroup/cpu.stat') or '').splitlines():
            k, _, v = line.partition(' ')
            if v.isdigit():
                cpu[k] = int(v)
        usage = cpu.get('usage_usec')
        ncpu = len(os.sched_getaffinity(0)) or 1        # CPUs this container may use (os.cpu_count() would report the whole host)
        pct = None
        if usage is not None and _last['usage'] is not None and now > _last['t']:
            pct = max(0.0, min(100.0, 100.0 * (usage - _last['usage']) / 1e6 / (now - _last['t']) / ncpu))
        try:
            free = shutil.disk_usage(config.JOBS_DIR if os.path.isdir(config.JOBS_DIR) else '/').free
        except OSError:
            free = None
        res = dict(mem_used_gb=round(stat.get('anon', 0) / 1e9, 1), mem_cache_gb=round(stat.get('file', 0) / 1e9, 1),
                   mem_limit_gb=round(int(limit) / 1e9, 1) if limit.isdigit() else None, cpu_pct=None if pct is None else round(pct), cpus=ncpu,
                   disk_free_gb=None if free is None else round(free / 1e9))
        _last.update(t=now, usage=usage, res=res)
        return res
