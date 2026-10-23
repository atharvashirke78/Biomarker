"""Server-side rendering of the live view (workflow chart, per-sample matrix, step table, resources) as HTML/SVG fragments.
Doing this in Python keeps the page's JavaScript to 'fetch and show' and makes the drawing testable."""
import html
import time

from . import resources, stages

# box positions of the chart (x, y), 160 x 40 each
POS = {'verify': (10, 10), 'fastp': (190, 10), 'star': (370, 10), 'sort': (550, 10), 'fc': (730, 10), 'salmon': (370, 70), 'qc': (10, 130),
       'merge': (190, 130), 'cohort': (370, 130), 'de': (550, 130), 'path': (730, 130), 'cand': (10, 190), 'model': (190, 190), 'report': (370, 190)}
EDGES = [('verify', 'fastp'), ('fastp', 'star'), ('fastp', 'salmon'), ('star', 'sort'), ('sort', 'fc'), ('fc', 'qc'), ('salmon', 'qc'), ('qc', 'merge'),
         ('merge', 'cohort'), ('cohort', 'de'), ('de', 'path'), ('path', 'cand'), ('cand', 'model'), ('model', 'report')]
W, H = 160, 40
MATRIX = [('FASTP', 'Trim'), ('STAR_ALIGN', 'Align'), ('SAMTOOLS_SORT_INDEX', 'Sort'), ('FEATURECOUNTS', 'Count'), ('SALMON_QUANT', 'Salmon'), ('SAMPLE_QC_FLAGS', 'QC')]
e = html.escape


def overall_percent(sts, job_status):
    if job_status == 'succeeded':
        return 100
    fr = [min(1.0, s['done'] / s['expected']) for s in sts if s['expected']]
    return int(round(100 * sum(fr) / len(fr))) if fr else 0


def chart_svg(sts):
    by = {s['key']: s for s in sts}
    out = ['<svg viewBox="0 0 900 250" width="100%" role="img" aria-label="Pipeline workflow with live status">',
           '<defs><marker id="ar" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#8a93a3"/></marker></defs>']
    for a, b in EDGES:
        (ax, ay), (bx, by_) = POS[a], POS[b]
        if by_ > ay:
            x1, y1, x2, y2 = ax + W / 2, ay + H, bx + (W / 2 if bx <= ax else 0), by_ + (0 if bx <= ax else H / 2)
        else:
            x1, y1, x2, y2 = ax + W, ay + H / 2, bx, by_ + H / 2
        flow = ' flow' if by[b]['status'] == 'running' else ''
        out.append('<line class="edge%s" x1="%g" y1="%g" x2="%g" y2="%g" marker-end="url(#ar)"/>' % (flow, x1, y1, x2, y2))
    for k, (x, y) in POS.items():
        s = by[k]
        if s['status'] == 'pending':
            sub = 'waiting'
        elif s['expected']:
            sub = '%d/%d' % (s['done'], s['expected'])
        else:
            sub = 'done' if s['status'] == 'done' else '%d done' % s['done']
        frac = min(1.0, s['done'] / s['expected']) if s['expected'] else (1.0 if s['status'] == 'done' else 0.0)
        out.append('<g class="st-%s"><title>%s: %s</title><rect class="box" x="%d" y="%d" width="%d" height="%d" rx="6"/>'
                   '<text x="%d" y="%d">%s</text><text class="sub" x="%d" y="%d">%s</text><rect class="bar" x="%d" y="%d" width="%.1f" height="3"/></g>'
                   % (s['status'], e(s['label']), e(s['what']), x, y, W, H, x + 8, y + 16, e(s['label']), x + 8, y + 32, e(sub), x + 2, y + H - 4, (W - 4) * frac))
    out.append('</svg>')
    return ''.join(out)


def matrix_html(by_sample, groups, qc_status):
    if not by_sample:
        return '<p class="small">Per-sample progress appears when the first tasks start.</p>'
    rows = ['<table class="matrix"><thead><tr><th>Sample</th><th>Group</th>' + ''.join('<th>%s</th>' % e(lbl) for _, lbl in MATRIX) + '<th>QC</th></tr></thead><tbody>']
    for sid in sorted(by_sample):
        st = by_sample[sid]
        cells = []
        for proc, lbl in MATRIX[:-1]:
            v = st.get(proc, 'pending')
            cells.append('<td class="c-%s" title="%s: %s">%s</td>' % (v, e(lbl), v, {'done': '&#10003;', 'running': '&#8226;&#8226;&#8226;', 'failed': '&#10007;'}.get(v, '')))
        q = qc_status.get(sid)
        cells.append('<td class="c-%s" title="QC flags">%s</td>' % ('done' if st.get('SAMPLE_QC_FLAGS') == 'done' else 'pending', '&#10003;' if st.get('SAMPLE_QC_FLAGS') == 'done' else ''))
        rows.append('<tr><td>%s</td><td>%s</td>%s<td class="q-%s">%s</td></tr>' % (e(sid), e(groups.get(sid, '')), ''.join(cells), e(q or ''), e(q or '')))
    return ''.join(rows) + '</tbody></table>'


def steps_html(sts):
    rows = ['<table class="steps"><thead><tr><th>Step</th><th>What it does</th><th>Tool</th><th>Status</th></tr></thead><tbody>']
    for s in sts:
        prog = ('%d/%d' % (s['done'], s['expected'])) if s['expected'] and s['status'] != 'pending' else ''
        rows.append('<tr class="r-%s"><td>%s</td><td>%s</td><td>%s</td><td><span class="pill p-%s">%s</span> %s</td></tr>'
                    % (s['status'], e(s['label']), e(s['what']), e(s['tool']), s['status'], s['status'], e(prog)))
    return ''.join(rows) + '</tbody></table>'


def resources_html(r):
    def bar(label, frac, text):
        f = 0 if frac is None else max(0.0, min(1.0, frac))
        return '<div class="gauge"><div class="gl">%s</div><div class="gb"><span style="width:%d%%"></span></div><div class="gt">%s</div></div>' % (e(label), round(100 * f), e(text))
    lim = r['mem_limit_gb']
    return ('<div class="gauges">' + bar('Memory in use by programs', (r['mem_used_gb'] / lim) if lim else None, '%.1f of %s GB' % (r['mem_used_gb'], lim if lim else '?'))
            + bar('CPU', None if r['cpu_pct'] is None else r['cpu_pct'] / 100, ('%d%% of %d CPUs' % (r['cpu_pct'], r['cpus'])) if r['cpu_pct'] is not None else '%d CPUs' % r['cpus'])
            + '<div class="small">File cache %.1f GB (reclaimable) · free disk %s GB</div></div>' % (r['mem_cache_gb'], r['disk_free_gb'] if r['disk_free_gb'] is not None else '?'))


def elapsed_text(started, finished, status):
    if not started:
        return 'not started'
    sec = int((finished or time.time()) - started)
    return '%d:%02d:%02d%s' % (sec // 3600, sec % 3600 // 60, sec % 60, '' if status != 'running' else ' (running)')


def live_html(job, view, groups=None, qc_status=None):
    """job: db row; view: server.job_view result (with stages and progress)."""
    sts = view.get('stages') or stages.compute({}, None, job['status'])
    prog = view.get('progress') or {}
    pct = overall_percent(sts, job['status'])
    running = [s for s in sts if s['status'] == 'running']
    if job['status'] == 'running' and running:
        now = '<b>Now:</b> ' + ' &middot; '.join('<b>%s</b> &ndash; %s' % (e(s['label']), e(s['what'])) for s in running[:3])
    elif job['status'] == 'succeeded':
        now = '<b>Finished.</b> Open the Summary, Figures and the report below.'
    elif job['status'] == 'failed':
        now = '<b>Stopped with an error:</b> ' + e(job.get('message') or '')
    elif job['status'] == 'queued':
        now = '<b>Waiting in the queue.</b> One job runs at a time.'
    elif job['status'] == 'cancelled':
        now = '<b>Cancelled.</b> You can resume it.'
    else:
        now = 'Not started yet.'
    head = ('<div class="liveHead"><div class="pct"><div class="pbar"><span style="width:%d%%"></span></div><b>%d%%</b> <span class="small">(approximate: share of finished tasks per step)</span></div>'
            '<div class="small">Elapsed %s · %d tasks finished</div></div><div class="now">%s</div>'
            % (pct, pct, e(elapsed_text(job.get('started'), job.get('finished'), job['status'])), prog.get('completed', 0), now))
    return (head + chart_svg(sts) + '<h3>Samples</h3>' + matrix_html(prog.get('by_sample') or {}, groups or {}, qc_status or {})
            + '<h3>Pod resources</h3>' + resources_html(resources.snapshot()) + '<h3>What each step does</h3>' + steps_html(sts))
