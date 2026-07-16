#!/usr/bin/env python3
"""Minimal web UI for driving telotron slurm runs and browsing results.

Runs a single-file stdlib HTTP server (no Flask). Use it from the login node
or an interactive slurm allocation:

    ./slurm/webui.py                 # binds 0.0.0.0:8765
    ./slurm/webui.py --port 9000

Endpoints:
  GET  /                    dashboard (squeue + DAG summary + results + rules)
  GET  /log/<name>          plain-text tail of the newest log matching <name>
  GET  /results/<path>      static-serve files under work/results/ (read-only)
  POST /run                 launch snakemake via slurm/submit.sh; form: target=<rule>
  POST /cancel              scancel every job named telo.* for this user

Nothing here bypasses slurm — every launch goes through `slurm/submit.sh`, so
per-rule cpus/mem/time come from `profiles/slurm/config.yaml`. Log tails are
truncated to the last 400 lines to keep responses bounded.
"""
from __future__ import annotations
import argparse
import glob
import html
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
os.chdir(REPO)
LOG_DIR = REPO / "work/logs/slurm"
RESULTS = REPO / "work/results"
WEB_LOG = REPO / "work/logs/webui.log"
LOG_DIR.mkdir(parents=True, exist_ok=True)
WEB_LOG.parent.mkdir(parents=True, exist_ok=True)

# ── Trusted rule set. Only these names can be POSTed to /run — prevents shell
#    injection via the form field. Extracted from Snakefile at import time.
def _discover_rules():
    txt = (REPO / "Snakefile").read_text()
    rules = []
    for line in txt.splitlines():
        line = line.strip()
        if line.startswith("rule ") and line.endswith(":"):
            rules.append(line[5:-1].strip())
        elif line.startswith("checkpoint ") and line.endswith(":"):
            rules.append(line[11:-1].strip())
    rules.append("all")
    return sorted(set(rules))

RULES = _discover_rules()
USER = os.environ.get("USER", "unknown")


def _run(cmd, timeout=15):
    """Run cmd (list), return (rc, stdout). Never raises."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode, (r.stdout or "") + (r.stderr or "")
    except FileNotFoundError:
        return 127, f"{cmd[0]}: not found\n"
    except subprocess.TimeoutExpired:
        return 124, f"{' '.join(cmd)}: timeout\n"


def squeue_html():
    if not shutil.which("squeue"):
        return "<p class='muted'>squeue unavailable on this host</p>"
    _, out = _run(["squeue", "-u", USER,
                   "-o", "%.10i %.30j %.10T %.10M %.6D %R"], timeout=10)
    if not out.strip():
        return "<p class='muted'>no jobs in queue</p>"
    return "<pre>" + html.escape(out.rstrip()) + "</pre>"


def dag_summary_html():
    snake = shutil.which("snakemake")
    if not snake:
        return "<p class='muted'>snakemake not on PATH — activate env before starting webui</p>"
    _, out = _run([snake, "--quiet", "-n", "--summary"], timeout=30)
    if not out.strip():
        return "<p class='muted'>(no DAG state yet)</p>"
    rows = out.splitlines()
    # Only show rows whose "status" col isn't ok — the interesting ones.
    keep = [rows[0]] if rows else []
    for r in rows[1:]:
        cols = r.split("\t")
        if len(cols) >= 4 and cols[3].strip() != "ok":
            keep.append(r)
    if len(keep) <= 1:
        return "<p class='muted'>every file up-to-date</p>"
    return "<pre>" + html.escape("\n".join(keep[:60])) + "</pre>"


def results_html():
    if not RESULTS.exists():
        return "<p class='muted'>no work/results/ yet</p>"
    interesting = [
        "telotron_pipeline_outputs.zip",
        "confident_species.tsv",
        "final_telotron_set_architecture.tsv",
        "final_species_summary.tsv",
        "all_species_raw_summary.tsv",
        "all_repeat_introns.tsv",
        "distance_to_end.tsv",
        "boundary_kmer_enrichment.tsv",
        "candidates_preview.html",
        "confident_report.html",
        "figures/telotron_ortholog_loci.pdf",
        "tert_deep_homology/confirmed_tert.tsv",
    ]
    rows = []
    for rel in interesting:
        p = RESULTS / rel
        if p.exists():
            size = p.stat().st_size
            mtime = time.strftime("%Y-%m-%d %H:%M", time.localtime(p.stat().st_mtime))
            rows.append((rel, size, mtime))
    if not rows:
        return "<p class='muted'>no results yet — run <code>all</code> or a specific rule</p>"
    body = "".join(
        f"<tr><td><a href='/results/{html.escape(rel)}'>{html.escape(rel)}</a></td>"
        f"<td class='num'>{fmt_size(sz)}</td><td class='muted'>{mt}</td></tr>"
        for rel, sz, mt in rows
    )
    return f"<table><thead><tr><th>Path</th><th>Size</th><th>Modified</th></tr></thead><tbody>{body}</tbody></table>"


def fmt_size(n):
    for unit in ("B", "K", "M", "G", "T"):
        if n < 1024:
            return f"{n:.1f}{unit}" if unit != "B" else f"{n}{unit}"
        n /= 1024
    return f"{n:.1f}P"


def logs_html():
    entries = sorted(LOG_DIR.glob("*.out"), key=lambda p: p.stat().st_mtime, reverse=True)[:30]
    if not entries:
        return "<p class='muted'>no logs under work/logs/slurm/ yet</p>"
    body = "".join(
        f"<tr><td><a href='/log/{html.escape(p.name)}'>{html.escape(p.name)}</a></td>"
        f"<td class='muted'>{time.strftime('%m-%d %H:%M', time.localtime(p.stat().st_mtime))}</td>"
        f"<td class='num'>{fmt_size(p.stat().st_size)}</td></tr>"
        for p in entries
    )
    return f"<table><thead><tr><th>Log</th><th>Modified</th><th>Size</th></tr></thead><tbody>{body}</tbody></table>"


def render_dashboard():
    rules_opts = "".join(f"<option value='{html.escape(r)}'>{html.escape(r)}</option>" for r in RULES)
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    return DASHBOARD_TEMPLATE.format(
        squeue=squeue_html(),
        dag=dag_summary_html(),
        results=results_html(),
        logs=logs_html(),
        rules=rules_opts,
        user=html.escape(USER),
        now=now,
    )


DASHBOARD_TEMPLATE = """\
<!doctype html>
<meta charset='utf-8'>
<meta name='viewport' content='width=device-width,initial-scale=1'>
<title>telotrons — pipeline UI</title>
<style>
  :root {{
    --bg:#0f1115; --card:#171a21; --border:#262a33; --fg:#e6e6e6;
    --muted:#8a8f99; --accent:#7aa2f7; --ok:#9ece6a; --warn:#e0af68;
    --err:#f7768e; --mono: ui-monospace,SFMono-Regular,Consolas,monospace;
  }}
  @media (prefers-color-scheme: light) {{
    :root {{ --bg:#f7f8fa; --card:#fff; --border:#e2e5eb; --fg:#1a1d24;
      --muted:#5e6470; --accent:#3b5bdb; --ok:#37864a; --warn:#a3671b; --err:#c92a2a; }}
  }}
  body {{ background:var(--bg); color:var(--fg); font:14px/1.5 system-ui,sans-serif;
    margin:0; padding:24px; max-width:1200px; margin:auto; }}
  h1 {{ font-size:20px; margin:0 0 4px; letter-spacing:-.01em; }}
  h2 {{ font-size:14px; margin:0 0 12px; color:var(--muted);
    text-transform:uppercase; letter-spacing:.06em; font-weight:600; }}
  header {{ display:flex; justify-content:space-between; align-items:baseline;
    margin-bottom:20px; padding-bottom:12px; border-bottom:1px solid var(--border); }}
  header .meta {{ color:var(--muted); font-family:var(--mono); font-size:12px; }}
  .grid {{ display:grid; gap:16px; grid-template-columns:repeat(auto-fit,minmax(360px,1fr)); }}
  .card {{ background:var(--card); border:1px solid var(--border); border-radius:8px;
    padding:16px; }}
  .card.full {{ grid-column:1/-1; }}
  pre {{ font-family:var(--mono); font-size:12px; background:var(--bg);
    border:1px solid var(--border); border-radius:6px; padding:10px;
    overflow:auto; max-height:280px; margin:0; }}
  table {{ width:100%; border-collapse:collapse; font-size:13px; }}
  th, td {{ text-align:left; padding:6px 8px; border-bottom:1px solid var(--border); }}
  th {{ font-weight:600; color:var(--muted); font-size:11px;
    text-transform:uppercase; letter-spacing:.05em; }}
  td.num {{ font-family:var(--mono); text-align:right; }}
  .muted {{ color:var(--muted); }}
  a {{ color:var(--accent); text-decoration:none; }}
  a:hover {{ text-decoration:underline; }}
  form {{ display:flex; gap:8px; align-items:center; flex-wrap:wrap; }}
  select, button {{ background:var(--card); color:var(--fg); border:1px solid var(--border);
    border-radius:6px; padding:6px 12px; font:inherit; cursor:pointer; }}
  button {{ background:var(--accent); color:#fff; border-color:var(--accent); font-weight:600; }}
  button.warn {{ background:var(--err); border-color:var(--err); }}
  button:hover {{ opacity:.9; }}
  .refresh {{ font-size:12px; }}
</style>

<header>
  <div>
    <h1>telotrons — pipeline UI</h1>
    <div class='meta'>user={user} · {now} · <a href='/' class='refresh'>refresh</a> ·
      auto-refresh in <span id='r'>30</span>s</div>
  </div>
  <form method='post' action='/run'>
    <select name='target'>{rules}</select>
    <button type='submit'>run</button>
  </form>
</header>

<div class='grid'>
  <div class='card'>
    <h2>slurm queue</h2>
    {squeue}
  </div>
  <div class='card'>
    <h2>DAG (needing update)</h2>
    {dag}
  </div>
  <div class='card full'>
    <h2>results</h2>
    {results}
  </div>
  <div class='card full'>
    <h2>recent logs</h2>
    {logs}
  </div>
  <div class='card full'>
    <h2>controls</h2>
    <form method='post' action='/cancel' onsubmit='return confirm("scancel every telo.* job?")'>
      <button class='warn' type='submit'>cancel all telo.* jobs</button>
    </form>
  </div>
</div>

<script>
  let s=30, el=document.getElementById('r');
  setInterval(()=>{{s--;el.textContent=s;if(s<=0)location.reload();}},1000);
</script>
"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        with open(WEB_LOG, "a") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {self.client_address[0]} {fmt % args}\n")

    def _send(self, code, body, ctype="text/html; charset=utf-8"):
        body_b = body.encode() if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body_b)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body_b)

    def _redirect(self, path="/"):
        self.send_response(303)
        self.send_header("Location", path)
        self.end_headers()

    def do_GET(self):
        path = urllib.parse.urlparse(self.path).path
        if path == "/":
            return self._send(200, render_dashboard())
        if path.startswith("/log/"):
            return self._serve_log(path[len("/log/"):])
        if path.startswith("/results/"):
            return self._serve_result(path[len("/results/"):])
        if path == "/health":
            return self._send(200, "ok\n", "text/plain")
        return self._send(404, "not found\n", "text/plain")

    def _serve_log(self, name):
        name = urllib.parse.unquote(name)
        # No path traversal — resolve, then check containment.
        p = (LOG_DIR / name).resolve()
        if LOG_DIR.resolve() not in p.parents and p != LOG_DIR.resolve():
            return self._send(400, "bad path\n", "text/plain")
        if not p.exists():
            return self._send(404, "no such log\n", "text/plain")
        # Tail last 400 lines to keep responses bounded.
        lines = p.read_text(errors="replace").splitlines()
        tail = "\n".join(lines[-400:])
        return self._send(200, tail, "text/plain; charset=utf-8")

    def _serve_result(self, rel):
        rel = urllib.parse.unquote(rel)
        p = (RESULTS / rel).resolve()
        if RESULTS.resolve() not in p.parents and p != RESULTS.resolve():
            return self._send(400, "bad path\n", "text/plain")
        if not p.exists() or not p.is_file():
            return self._send(404, "not found\n", "text/plain")
        # Content-type sniff by suffix only — no magic bytes needed for this use.
        ext = p.suffix.lower()
        ctype = {
            ".html": "text/html; charset=utf-8",
            ".tsv":  "text/tab-separated-values",
            ".csv":  "text/csv",
            ".txt":  "text/plain; charset=utf-8",
            ".pdf":  "application/pdf",
            ".png":  "image/png",
            ".svg":  "image/svg+xml",
            ".zip":  "application/zip",
            ".json": "application/json",
        }.get(ext, "application/octet-stream")
        # Files can be large (multi-GB zip); stream by reading in chunks.
        size = p.stat().st_size
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(size))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        with open(p, "rb") as fh:
            while True:
                chunk = fh.read(65536)
                if not chunk:
                    break
                try:
                    self.wfile.write(chunk)
                except BrokenPipeError:
                    return

    def do_POST(self):
        length = int(self.headers.get("Content-Length", "0") or 0)
        raw = self.rfile.read(length).decode("utf-8", errors="replace")
        form = urllib.parse.parse_qs(raw)
        path = urllib.parse.urlparse(self.path).path
        if path == "/run":
            return self._post_run(form)
        if path == "/cancel":
            return self._post_cancel()
        return self._send(404, "not found\n", "text/plain")

    def _post_run(self, form):
        target = (form.get("target") or [""])[0].strip()
        # Whitelist: only rules discovered from the Snakefile can be launched.
        if target not in RULES:
            return self._send(400, f"unknown target: {target!r}\n", "text/plain")
        # Detach: run submit.sh in the background so the driver survives the
        # HTTP round-trip. Each launch gets its own timestamped log so multi
        # driver instances don't stomp each other.
        stamp = time.strftime("%Y%m%d-%H%M%S")
        out = LOG_DIR / f"driver.{target}.{stamp}.out"
        cmd = f"nohup ./slurm/submit.sh {shlex.quote(target)} > {shlex.quote(str(out))} 2>&1 &"
        subprocess.Popen(["bash", "-lc", cmd], cwd=str(REPO), start_new_session=True)
        with open(WEB_LOG, "a") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} LAUNCH {target} log={out.name}\n")
        return self._redirect("/")

    def _post_cancel(self):
        if not shutil.which("scancel"):
            return self._send(500, "scancel not available\n", "text/plain")
        # Cancel by job-name pattern rather than a scancel-all.
        rc, out = _run(["scancel", "-u", USER, "--name=telo.*"], timeout=15)
        with open(WEB_LOG, "a") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} CANCEL rc={rc}\n")
        return self._redirect("/")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--port", type=int, default=8765)
    args = ap.parse_args()
    srv = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"telotrons webui: http://{args.host}:{args.port}  (user={USER}, repo={REPO})",
          file=sys.stderr, flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        srv.shutdown()


if __name__ == "__main__":
    main()
