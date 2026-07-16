#!/usr/bin/env python3
"""Dynamic report tracing the confident telotron-bearing species set.

Sections:
  1. Provenance funnel — pre-filter → filter → dedup → confident.
  2. Confident-species table — per-species stats + annotation slots for
     downstream analyses (TERT +/-, expression, ortholog fill/birth, ...).
  3. Excluded (admitted but below confident bar) — small section.
  4. Per-species detail — telotrons with 50-bp flanks, motif-direction
     coloring, and architecture badges.

Update semantics: re-run to overwrite work/results/confident_report.html.
As downstream stages complete, pass their outputs via CLI flags and they
appear as annotations in the species table without changing the URL.
"""
import argparse
import csv
import gzip
import html
import json
import os
import sys
from collections import defaultdict, Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import find_genome_fasta  # noqa: E402

FLANK = 50
RC_TAB = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(RC_TAB)[::-1]


def rotations(m: str) -> set:
    return {m[i:] + m[:i] for i in range(len(m))}


def load_fasta(path: str) -> dict:
    opener = gzip.open if path.endswith(".gz") else open
    seqs, cur, buf = {}, None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf).upper()
    return seqs


def build_kmer_sets(motif: str):
    L = len(motif)
    fwd = rotations(motif)
    rev = rotations(rc(motif))

    def mismatches(kms):
        out = set()
        for km in kms:
            for i in range(L):
                for b in "ACGT":
                    if b != km[i]:
                        out.add(km[:i] + b + km[i + 1:])
        return out

    return fwd, mismatches(fwd) - fwd - rev, rev, mismatches(rev) - rev - fwd


CLS_RANK = {"ge": 5, "ce": 4, "gm": 3, "cm": 2}


def classify_positions(seq: str, motif: str):
    n = len(seq); L = len(motif)
    if L == 0 or n < L:
        return [None] * n
    fwd_e, fwd_m, rev_e, rev_m = build_kmer_sets(motif)
    cls = [None] * n
    for i in range(n - L + 1):
        km = seq[i:i + L]
        if km in fwd_e:
            new = "ge"
        elif km in rev_e:
            new = "ce"
        elif km in fwd_m:
            new = "gm"
        elif km in rev_m:
            new = "cm"
        else:
            continue
        for j in range(i, i + L):
            if cls[j] is None or CLS_RANK[new] > CLS_RANK[cls[j]]:
                cls[j] = new
    return cls


def render_seq(seq: str, motif: str, up_len: int, dn_len: int) -> str:
    n = len(seq)
    intron_start = up_len
    intron_end = n - dn_len
    cls = classify_positions(seq, motif)
    eff = []
    for i in range(n):
        c = cls[i] or ("fl" if (i < intron_start or i >= intron_end) else "in")
        eff.append(c)
    parts = []
    buf, cur = [], None
    for i in range(n):
        if i == intron_start:
            if buf: parts.append(f'<span class="{cur}">{"".join(buf)}</span>'); buf, cur = [], None
            parts.append('<span class="sep">⎢</span>')
        if i == intron_end and intron_end > intron_start:
            if buf: parts.append(f'<span class="{cur}">{"".join(buf)}</span>'); buf, cur = [], None
            parts.append('<span class="sep">⎢</span>')
        c = eff[i]
        if c != cur:
            if buf: parts.append(f'<span class="{cur}">{"".join(buf)}</span>'); buf = []
            cur = c
        buf.append(html.escape(seq[i]))
    if buf:
        parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
    return "".join(parts)


CSS = """
:root {
  --bg: #fbfaf7; --fg: #16161a; --muted: #6a6a72;
  --card: #ffffff; --line: #e6e4dd; --accent: #1b5f3a;
  --g-exact: #1b7a3a; --g-mm: #7fb992;
  --c-exact: #b2182b; --c-mm: #d68b9a;
  --flank: #999; --intron: #2b2b2f; --sep: #a89f7c;
  --chip-bg: #efece5;
  --bar-bg: #eae7de; --bar-fg: #6b9b7f;
  --arch-fr: #1b5f3a; --arch-f: #2f78a8; --arch-r: #a05a2c;
  --arch-linker: #8b4a8a; --arch-multi: #7a7a4a; --arch-other: #7b7b83;
  --status-pass: #1b7a3a; --status-fail: #b2182b; --status-pending: #8a8a91;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #17181b; --fg: #ececee; --muted: #8a8a91;
    --card: #1e2024; --line: #2a2c31; --accent: #7ec89d;
    --g-exact: #5ec97e; --g-mm: #a9d5b9;
    --c-exact: #ff7683; --c-mm: #f0b6bd;
    --flank: #7d7d85; --intron: #c9c9cd; --sep: #c2b477;
    --chip-bg: #262930;
    --bar-bg: #2a2c31; --bar-fg: #6bb389;
    --arch-fr: #5ec97e; --arch-f: #74b0dc; --arch-r: #d99866;
    --arch-linker: #c489bf; --arch-multi: #bcbc82; --arch-other: #9ea0a8;
    --status-pass: #5ec97e; --status-fail: #ff7683; --status-pending: #7a7a80;
  }
}
:root[data-theme="light"] {
  --bg: #fbfaf7; --fg: #16161a; --muted: #6a6a72;
  --card: #ffffff; --line: #e6e4dd; --accent: #1b5f3a;
  --g-exact: #1b7a3a; --g-mm: #7fb992;
  --c-exact: #b2182b; --c-mm: #d68b9a;
  --flank: #999; --intron: #2b2b2f; --sep: #a89f7c;
  --chip-bg: #efece5;
  --bar-bg: #eae7de; --bar-fg: #6b9b7f;
  --arch-fr: #1b5f3a; --arch-f: #2f78a8; --arch-r: #a05a2c;
  --arch-linker: #8b4a8a; --arch-multi: #7a7a4a; --arch-other: #7b7b83;
  --status-pass: #1b7a3a; --status-fail: #b2182b; --status-pending: #8a8a91;
}
:root[data-theme="dark"] {
  --bg: #17181b; --fg: #ececee; --muted: #8a8a91;
  --card: #1e2024; --line: #2a2c31; --accent: #7ec89d;
  --g-exact: #5ec97e; --g-mm: #a9d5b9;
  --c-exact: #ff7683; --c-mm: #f0b6bd;
  --flank: #7d7d85; --intron: #c9c9cd; --sep: #c2b477;
  --chip-bg: #262930;
  --bar-bg: #2a2c31; --bar-fg: #6bb389;
  --arch-fr: #5ec97e; --arch-f: #74b0dc; --arch-r: #d99866;
  --arch-linker: #c489bf; --arch-multi: #bcbc82; --arch-other: #9ea0a8;
  --status-pass: #5ec97e; --status-fail: #ff7683; --status-pending: #7a7a80;
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--fg);
  font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
  font-size: 14px; line-height: 1.5;
}
.container { max-width: 1400px; margin: 0 auto; padding: 0 32px; }
header.top {
  padding: 28px 32px 20px; background: var(--card);
  border-bottom: 1px solid var(--line);
}
h1 { font-size: 22px; font-weight: 600; margin: 0 0 4px; letter-spacing: -0.01em; }
h2 { font-size: 16px; font-weight: 600; margin: 32px 0 10px; letter-spacing: -0.005em; }
h3 { font-size: 14px; font-weight: 600; margin: 20px 0 4px;
     display: flex; align-items: baseline; gap: 12px; padding-top: 8px; }
h3 .org { font-weight: 400; font-style: italic; color: var(--muted); font-size: 12.5px; }
p.meta { color: var(--muted); font-size: 12.5px; margin: 4px 0 8px; }
.subtitle { font-size: 13px; color: var(--muted); margin: 0 0 12px; }
/* Funnel */
.funnel {
  display: flex; align-items: stretch; gap: 4px; margin: 16px 0 4px; flex-wrap: wrap;
}
.funnel .step {
  padding: 10px 14px; border: 1px solid var(--line); border-radius: 6px;
  background: var(--bg); min-width: 140px;
}
.funnel .step .n {
  font-size: 20px; font-weight: 600; letter-spacing: -0.01em;
  font-variant-numeric: tabular-nums;
}
.funnel .step .lbl { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.06em; }
.funnel .step .sub { font-size: 11.5px; color: var(--muted); margin-top: 2px; }
.funnel .arrow { display: flex; align-items: center; color: var(--muted); font-size: 18px; padding: 0 4px; }
.funnel .step.confident { border-color: var(--accent); background: var(--card); }
.funnel .step.confident .n { color: var(--accent); }
/* Legend */
.legend {
  display: flex; flex-wrap: wrap; gap: 12px 20px; align-items: center;
  margin: 12px 0 4px; font-size: 11.5px; color: var(--muted);
  font-family: ui-monospace, monospace;
}
.legend .sw { display: inline-block; width: 12px; height: 12px;
              vertical-align: -2px; margin-right: 4px; border-radius: 2px; }
/* Species table */
table.species {
  width: 100%; border-collapse: collapse; margin: 8px 0 24px;
  font-size: 12.5px; font-variant-numeric: tabular-nums;
}
table.species th {
  text-align: left; padding: 8px 10px 6px; font-weight: 600;
  color: var(--muted); font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em;
  border-bottom: 1px solid var(--line);
}
table.species th.numeric, table.species td.numeric { text-align: right; }
table.species tr { border-bottom: 1px solid var(--line); }
table.species tr:hover { background: var(--chip-bg); }
table.species td {
  padding: 8px 10px; vertical-align: middle;
}
table.species td.gid { font-family: ui-monospace, monospace; font-weight: 500; }
table.species td.gid a { color: var(--fg); text-decoration: none; }
table.species td.gid a:hover { color: var(--accent); }
table.species td.org { font-style: italic; color: var(--muted); }
.bar { display: inline-block; width: 60px; height: 10px; background: var(--bar-bg);
       border-radius: 2px; vertical-align: middle; overflow: hidden; }
.bar > i {
  display: block; height: 100%; background: var(--bar-fg); border-radius: 2px;
}
.status { font-family: ui-monospace, monospace; font-size: 11px; }
.status.pass { color: var(--status-pass); }
.status.fail { color: var(--status-fail); }
.status.pending { color: var(--status-pending); }
.archbar {
  display: inline-flex; height: 10px; border-radius: 2px; overflow: hidden;
  min-width: 60px; vertical-align: middle;
}
.archbar > span { height: 100%; display: inline-block; }
.archbar > span.fr { background: var(--arch-fr); }
.archbar > span.f { background: var(--arch-f); }
.archbar > span.r { background: var(--arch-r); }
.archbar > span.l { background: var(--arch-linker); }
.archbar > span.m { background: var(--arch-multi); }
.archbar > span.o { background: var(--arch-other); }
.arch-badge {
  display: inline-block; padding: 1px 5px; border-radius: 2px;
  font-size: 10px; font-family: ui-monospace, monospace; text-transform: uppercase;
  letter-spacing: 0.03em; color: var(--card); font-weight: 500;
}
.arch-badge.fr { background: var(--arch-fr); }
.arch-badge.f { background: var(--arch-f); }
.arch-badge.r { background: var(--arch-r); }
.arch-badge.l { background: var(--arch-linker); }
.arch-badge.m { background: var(--arch-multi); }
.arch-badge.o { background: var(--arch-other); }
/* Locus rows */
.row { padding: 8px 0; border-bottom: 1px solid var(--line); }
.rmeta {
  font-size: 11.5px; display: flex; flex-wrap: wrap;
  gap: 4px 10px; align-items: baseline;
}
.lbl { font-family: ui-monospace, monospace; color: var(--muted); }
.chip {
  background: var(--chip-bg); padding: 1px 6px; border-radius: 3px;
  font-size: 10.5px; color: var(--muted); font-variant-numeric: tabular-nums;
}
.chip.pur { color: var(--g-exact); font-weight: 600; }
pre.seq {
  margin: 4px 0 0; padding: 6px 10px;
  background: var(--bg); border: 1px solid var(--line); border-radius: 4px;
  font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
  font-size: 12.5px; letter-spacing: 0.03em;
  overflow-x: auto; white-space: pre;
}
.seq .ge { color: var(--g-exact); font-weight: 600; }
.seq .gm { color: var(--g-mm); }
.seq .ce { color: var(--c-exact); font-weight: 600; }
.seq .cm { color: var(--c-mm); }
.seq .fl { color: var(--flank); }
.seq .in { color: var(--intron); }
.seq .sep { color: var(--sep); margin: 0 1px; font-weight: 600; }
.missing { color: var(--muted); font-style: italic; padding: 6px 0; }
.excluded {
  margin-top: 20px; padding: 12px 16px; background: var(--chip-bg); border-radius: 6px;
  font-size: 12.5px; color: var(--muted);
}
.excluded strong { color: var(--fg); font-weight: 600; }
details.filters {
  margin: 16px 0 8px; padding: 8px 14px;
  background: var(--card); border: 1px solid var(--line); border-radius: 6px;
}
details.filters summary {
  font-weight: 600; font-size: 13px; cursor: pointer;
  padding: 4px 0; color: var(--fg); list-style: none;
}
details.filters summary::before { content: "▸ "; color: var(--muted); }
details.filters[open] summary::before { content: "▾ "; color: var(--muted); }
.filter-grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 12px; margin-top: 10px;
}
.fstage {
  padding: 10px 14px; background: var(--bg); border: 1px solid var(--line); border-radius: 4px;
  font-size: 12px;
}
.fstage h4 { margin: 0 0 6px; font-size: 12.5px; font-weight: 600; }
.fstage h4 .muted { font-weight: 400; color: var(--muted); font-size: 11.5px; }
.fstage ul { margin: 4px 0 0; padding-left: 16px; }
.fstage li { margin: 2px 0; line-height: 1.4; }
.fstage code {
  background: var(--chip-bg); padding: 1px 4px; border-radius: 2px;
  font-size: 11px; font-family: ui-monospace, monospace;
}
.fstage em { color: var(--muted); font-style: italic; }
"""


ARCH_KEY = {
    "GT-F-R-AG": ("fr", "F-R"),
    "GT-F-AG": ("f", "F"),
    "GT-R-AG": ("r", "R"),
    "GT-R-linker-F-AG": ("l", "R-L-F"),
    "GT-F-linker-R-AG": ("l", "F-L-R"),
    "Multi-junction": ("m", "MULTI"),
    "Other": ("o", "OTHER"),
}


def arch_key(a): return ARCH_KEY.get(a, ("o", a or "?"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="work/results/final_telotron_set_architecture.tsv")
    ap.add_argument("--confident", default="work/results/confident_species.tsv")
    ap.add_argument("--pre-filter-summary", default="work/results/all_species_raw_summary.tsv")
    ap.add_argument("--filter-final", default="work/results/final_telotron_set.tsv")
    ap.add_argument("--dedup", default="work/results/final_telotron_set_dedup.tsv")
    ap.add_argument("--canonical", default="work/manifests/canonical_motifs.tsv")
    ap.add_argument("--refseq-dir", default="data/raw/refseq")
    ap.add_argument("--tara-dir", default="data/raw/tara")
    ap.add_argument("--out", default="work/results/confident_report.html")
    ap.add_argument("--max-per-species", type=int, default=60)
    # Annotation slots (JSON files keyed by genome_id) — passed as
    # `--annotate NAME=path/to/data.json`.
    ap.add_argument("--annotate", action="append", default=[],
                    help="NAME=json path; JSON maps genome_id -> string label")
    args = ap.parse_args()

    # ── provenance ──────────────────────────────────────────────
    pre_cands = 0; pre_species = 0
    if os.path.exists(args.pre_filter_summary):
        for r in csv.DictReader(open(args.pre_filter_summary), delimiter="\t"):
            pre_species += 1
            try: pre_cands += int(r.get("telotron_candidates", 0) or 0)
            except: pass
    n_filter = sum(1 for _ in open(args.filter_final)) - 1 if os.path.exists(args.filter_final) else 0
    filter_species = len({r["genome_id"] for r in csv.DictReader(open(args.filter_final), delimiter="\t")}) if os.path.exists(args.filter_final) else 0
    n_dedup = sum(1 for _ in open(args.dedup)) - 1 if os.path.exists(args.dedup) else n_filter
    dedup_species = len({r["genome_id"] for r in csv.DictReader(open(args.dedup), delimiter="\t")}) if os.path.exists(args.dedup) else filter_species

    # ── confident species ─────────────────────────────────────────
    confident = list(csv.DictReader(open(args.confident), delimiter="\t"))
    confident_ids = {r["genome_id"] for r in confident}
    n_conf_loci = sum(int(r["n_telotrons"]) for r in confident)

    # ── arch rows keyed by genome ─────────────────────────────────
    arch_rows_by_gid = defaultdict(list)
    for r in csv.DictReader(open(args.arch), delimiter="\t"):
        arch_rows_by_gid[r["genome_id"]].append(r)

    # ── canonical motif per species ───────────────────────────────
    canonical = {}
    if os.path.exists(args.canonical):
        for r in csv.DictReader(open(args.canonical), delimiter="\t"):
            if r.get("motif"): canonical[r["genome_id"]] = r["motif"]

    # ── annotations (per genome_id) ───────────────────────────────
    annotations = {}
    # Auto-annotations from standard paths (fills in as downstream stages complete).
    def _auto_tert():
        p = "work/results/tert_deep_homology/confirmed_tert.tsv"
        if not os.path.exists(p): return None
        counts = defaultdict(int); best_ev = {}
        for r in csv.DictReader(open(p), delimiter="\t"):
            gid = r.get("genome_id", "")
            counts[gid] += 1
            try:
                ev = float(r.get("min_ievalue", r.get("i_evalue", "1e30")))
                if gid not in best_ev or ev < best_ev[gid]:
                    best_ev[gid] = ev
            except (ValueError, TypeError):
                pass
        out = {}
        for gid, n in counts.items():
            ev = best_ev.get(gid)
            out[gid] = f"+ {n} (e={ev:.0e})" if ev is not None else f"+ {n}"
        return out
    def _auto_orthologs():
        # Summarize telotron_orthologs by counting FILL / BIRTH calls per focal genome.
        p = "work/results/telotron_orthologs/fill_vs_birth.tsv"
        if not os.path.exists(p): return None
        counts = defaultdict(lambda: {"FILL": 0, "BIRTH": 0, "ABS": 0, "?": 0})
        for r in csv.DictReader(open(p), delimiter="\t"):
            gid = r.get("focal_genome_id") or r.get("genome_id") or ""
            call = r.get("call", "?")
            counts[gid][call] = counts[gid].get(call, 0) + 1
        return {gid: f"F{c['FILL']}/B{c['BIRTH']}" for gid, c in counts.items()}
    def _auto_interstitial():
        p = "work/results/interstitial_arrays.tsv"
        if not os.path.exists(p): return None
        counts = defaultdict(int)
        for r in csv.DictReader(open(p), delimiter="\t"):
            counts[r.get("genome_id", "")] += 1
        return {gid: str(n) for gid, n in counts.items()}

    for name, fn in [("TERT", _auto_tert), ("Ortho F/B", _auto_orthologs), ("Interstitial", _auto_interstitial)]:
        got = fn()
        if got is not None:
            annotations[name] = got

    for spec in args.annotate:
        if "=" not in spec: continue
        name, path = spec.split("=", 1)
        try:
            annotations[name] = json.load(open(path))
        except Exception as e:
            print(f"[warn] annotation {name} @{path}: {e}", flush=True)
            annotations[name] = {}
    annot_names = list(annotations)

    # ── per-species architecture breakdown ────────────────────────
    def arch_breakdown(gid):
        c = Counter()
        for r in arch_rows_by_gid.get(gid, []):
            c[r.get("architecture", "?")] += 1
        return c

    def arch_bar_html(c, total):
        # order archs consistently
        order = ["GT-F-R-AG", "GT-F-AG", "GT-R-AG", "GT-R-linker-F-AG",
                 "GT-F-linker-R-AG", "Multi-junction", "Other"]
        segments = []
        for a in order:
            n = c.get(a, 0)
            if not n: continue
            k, _ = arch_key(a)
            pct = 100 * n / total
            segments.append(f'<span class="{k}" style="width:{pct:.1f}%" '
                            f'title="{a}: {n} ({pct:.0f}%)"></span>')
        return f'<div class="archbar">{"".join(segments)}</div>'

    # sort confident by n_telotrons desc
    confident.sort(key=lambda r: -int(r["n_telotrons"]))

    # ── build FASTA cache lazily ──────────────────────────────────
    fasta_cache = {}

    def get_seq(gid, seqid):
        if gid not in fasta_cache:
            p = find_genome_fasta(gid, args.refseq_dir, args.tara_dir, required=False)
            fasta_cache[gid] = load_fasta(p) if p else {}
        return fasta_cache[gid].get(seqid)

    # ── excluded (admitted but not confident) — from arch rows ─────
    excluded_by_gid = {gid: rows for gid, rows in arch_rows_by_gid.items() if gid not in confident_ids}

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        f.write('<meta charset="utf-8"><title>Telotron confident-set report</title>'
                f'<style>{CSS}</style></head><body>')
        f.write('<header class="top"><div class="container">')
        f.write('<h1>Telotron confident-set report</h1>')
        f.write('<p class="subtitle">Live-updating trace of the candidate set as it\'s pared down and annotated.</p>')

        # Funnel
        max_n = max(pre_cands, n_filter, n_dedup, n_conf_loci, 1)
        f.write('<div class="funnel">')
        f.write(f'<div class="step"><div class="lbl">Pre-filter scan</div>'
                f'<div class="n">{pre_cands:,}</div><div class="sub">{pre_species} species</div></div>')
        f.write('<div class="arrow">→</div>')
        f.write(f'<div class="step"><div class="lbl">filter_final</div>'
                f'<div class="n">{n_filter:,}</div><div class="sub">{filter_species} species '
                f'· motif gate OFF · single≥0.85 or bidir≥0.40+3</div></div>')
        f.write('<div class="arrow">→</div>')
        f.write(f'<div class="step"><div class="lbl">dedup</div>'
                f'<div class="n">{n_dedup:,}</div><div class="sub">{dedup_species} species '
                f'· both-flank blastn</div></div>')
        f.write('<div class="arrow">→</div>')
        f.write(f'<div class="step confident"><div class="lbl">Confident bearers</div>'
                f'<div class="n">{n_conf_loci:,}</div><div class="sub">{len(confident)} species '
                f'· n≥3 OR bidir≥2</div></div>')
        f.write('</div>')

        # Filters-applied panel
        f.write('<details class="filters" open><summary>Filters applied at each stage</summary>')
        f.write('<div class="filter-grid">')
        f.write('<div class="fstage"><h4>Stage 4 · scan_all <span class="muted">(per-intron sieve)</span></h4>'
                '<ul>'
                '<li><code>min_repeat_frac ≥ 0.30</code> — intron ≥30% telomeric-motif-covered</li>'
                '<li><code>max_flank_repeat_frac ≤ 0.25</code> — flanks NOT themselves telomeric (kills misannotated subtelomere introns)</li>'
                '<li><code>min_intron_len ≥ 30 bp</code> — kills 1-bp degenerate GFF introns</li>'
                '</ul></div>')
        f.write('<div class="fstage"><h4>Stage 5 · filter_final <span class="muted">(admission)</span></h4>'
                '<ul>'
                '<li><b>Single-array</b>: <code>telomeric_frac ≥ 0.85</code></li>'
                '<li><b>OR bidirectional</b>: <code>frac ≥ 0.40 AND fwd_hits ≥ 3 AND rev_hits ≥ 3</code></li>'
                '<li>+ <code>collapse_unique_loci=true</code> (alt-splice paralogs collapse)</li>'
                '<li>+ <code>require_terminal_motif_match=false</code> <em>(off this run)</em></li>'
                '<li>+ <code>require_canonical_splice=false</code> (GT-AG and non-GT-AG both kept)</li>'
                '</ul></div>')
        f.write('<div class="fstage"><h4>Stage 6 · dedup_telotrons <span class="muted">(cross-locus)</span></h4>'
                '<ul>'
                '<li>Within-species all-vs-all blastn of 250bp upstream + 250bp downstream flanks</li>'
                '<li>Duplicate iff <b>both</b> flanks share <code>≥100 bp @ ≥95% identity</code></li>'
                '<li>Longest-intron representative kept per cluster</li>'
                '</ul></div>')
        f.write('<div class="fstage"><h4>Stage 7 · confident_species <span class="muted">(species bar)</span></h4>'
                '<ul>'
                '<li><b>n_telotrons ≥ 3</b> (<code>min_n=3</code>)</li>'
                '<li><b>OR n_bidirectional ≥ 2</b> (<code>min_bidir=2</code>)</li>'
                '<li>Bidirectional architectures: <code>GT-F-R-AG</code>, <code>GT-R-linker-F-AG</code>, <code>GT-F-linker-R-AG</code>, <code>Multi-junction</code></li>'
                '</ul></div>')
        f.write('</div></details>')

        # Legend
        f.write('<div class="legend">'
                '<span><span class="sw" style="background:var(--g-exact)"></span>G-rich exact</span>'
                '<span><span class="sw" style="background:var(--g-mm)"></span>G-rich 1-mm</span>'
                '<span><span class="sw" style="background:var(--c-exact)"></span>C-rich exact</span>'
                '<span><span class="sw" style="background:var(--c-mm)"></span>C-rich 1-mm</span>'
                '<span><span class="sw" style="background:var(--flank)"></span>flank</span>'
                '<span><span class="sw" style="background:var(--intron)"></span>intron non-motif</span>'
                '<span style="border-left:1px solid var(--line);padding-left:12px">'
                'Arch: '
                '<span class="arch-badge fr">F-R</span> '
                '<span class="arch-badge f">F</span> '
                '<span class="arch-badge r">R</span> '
                '<span class="arch-badge l">LINK</span> '
                '<span class="arch-badge m">MULTI</span> '
                '<span class="arch-badge o">OTHER</span></span>'
                '</div>')
        f.write('</div></header>')

        f.write('<main class="container">')

        # Species table
        f.write(f'<h2>Confident species ({len(confident)}, {n_conf_loci:,} loci)</h2>')
        f.write('<table class="species"><thead><tr>')
        f.write('<th>Species</th><th>Organism</th><th>Source</th>')
        f.write('<th class="numeric">Telotrons</th>')
        f.write('<th class="numeric">Bidir</th>')
        f.write('<th>Architecture mix</th>')
        f.write('<th class="numeric">Med. telo_frac</th>')
        f.write('<th class="numeric">Med. len</th>')
        for name in annot_names:
            f.write(f'<th>{html.escape(name)}</th>')
        f.write('</tr></thead><tbody>')
        for r in confident:
            gid = r["genome_id"]
            slug = gid.replace(".", "_").replace("/", "_")
            n = int(r["n_telotrons"])
            nb = int(r["n_bidirectional"])
            arch_rows = arch_rows_by_gid.get(gid, [])
            bd = arch_breakdown(gid)
            tfracs = sorted(float(x.get("telomeric_frac", 0) or 0) for x in arch_rows)
            ilens = sorted(int(x.get("intron_len", 0) or 0) for x in arch_rows)
            med_tf = tfracs[len(tfracs) // 2] if tfracs else 0
            med_il = ilens[len(ilens) // 2] if ilens else 0
            bar_w = int(60 * n / max(n_conf_loci, 1) * len(confident))  # scale within table
            bar_w = min(60, bar_w)
            f.write(f'<tr><td class="gid"><a href="#s-{slug}">{html.escape(gid)}</a></td>'
                    f'<td class="org">{html.escape(r["organism"])}</td>'
                    f'<td>{html.escape(r["source"])}</td>'
                    f'<td class="numeric">{n:,} <span class="bar"><i style="width:{bar_w/60*100:.0f}%"></i></span></td>'
                    f'<td class="numeric">{nb:,}</td>'
                    f'<td>{arch_bar_html(bd, sum(bd.values()) or 1)}</td>'
                    f'<td class="numeric">{med_tf:.2f}</td>'
                    f'<td class="numeric">{med_il}</td>')
            for name in annot_names:
                v = annotations[name].get(gid, "")
                if v == "":
                    f.write('<td><span class="status pending">—</span></td>')
                else:
                    cls = "pass" if str(v).startswith("+") else ("fail" if str(v).startswith("-") else "pending")
                    f.write(f'<td><span class="status {cls}">{html.escape(str(v))}</span></td>')
            f.write('</tr>')
        f.write('</tbody></table>')

        # Excluded
        if excluded_by_gid:
            excl_counts = [(gid, len(rows)) for gid, rows in excluded_by_gid.items()]
            excl_counts.sort(key=lambda x: -x[1])
            f.write('<div class="excluded"><strong>Below confident bar '
                    f'({len(excl_counts)} species, {sum(c for _,c in excl_counts)} loci):</strong> ')
            items = []
            for gid, n in excl_counts:
                org = (arch_rows_by_gid[gid][0].get("organism") or "").split()[:2]
                items.append(f'<span class="lbl">{html.escape(gid)}</span> '
                             f'<span class="chip">n={n}</span> '
                             f'<span class="chip">{html.escape(" ".join(org))}</span>')
            f.write(" &nbsp; ".join(items))
            f.write('</div>')

        # Per-species detail
        f.write('<h2>Per-species detail</h2>')
        f.write(f'<p class="meta">Top {args.max_per_species} loci per species, sorted by '
                f'telomeric_frac desc → telomeric_bases desc → intron_len desc. '
                f'±{FLANK} bp flanks; architecture badge on each row.</p>')

        for r in confident:
            gid = r["genome_id"]
            slug = gid.replace(".", "_").replace("/", "_")
            org = r["organism"]
            arch_rows = list(arch_rows_by_gid.get(gid, []))
            arch_rows.sort(key=lambda x: (
                -float(x.get("telomeric_frac", 0) or 0),
                -int(float(x.get("telomeric_bases", 0) or 0)),
                -int(x.get("intron_len", 0) or 0),
            ))
            n = len(arch_rows); display = arch_rows[:args.max_per_species]
            capped = n > args.max_per_species
            motif_c = canonical.get(gid, "") or ""

            f.write(f'<h3 id="s-{slug}">{html.escape(gid)} '
                    f'<span class="org">{html.escape(org)}</span></h3>')
            note = f' (top {args.max_per_species} of {n})' if capped else ''
            canon = html.escape(motif_c) if motif_c else '<em>(none — per-locus motif)</em>'
            f.write(f'<p class="meta">{r["source"]} · {n} loci{note} · canonical motif: <code>{canon}</code></p>')

            for row in display:
                try:
                    start, end = int(row["start"]), int(row["end"])
                except (KeyError, ValueError):
                    continue
                seqid = row["seqid"]
                full = get_seq(gid, seqid)
                if not full:
                    f.write(f'<div class="row missing">{html.escape(seqid)}:{start}-{end} — seq not found</div>')
                    continue
                # GFF is 1-based inclusive; python slicing is 0-based half-open.
                # Intron occupies python[start-1 : end].
                pystart, pyend = start - 1, end
                up_s = max(0, pystart - FLANK); dn_e = min(len(full), pyend + FLANK)
                seg = full[up_s:dn_e]
                up_len = pystart - up_s; dn_len = dn_e - pyend
                strand = row.get("strand", "+")
                if strand == "-":
                    seg = rc(seg); up_len, dn_len = dn_len, up_len
                motif = (row.get("motif") or motif_c or "TTAGGG").strip()
                colored = render_seq(seg, motif, up_len, dn_len)
                try: tfrac = float(row["telomeric_frac"])
                except: tfrac = 0.0
                try: tbases = int(float(row.get("telomeric_bases", 0)))
                except: tbases = 0
                try: ilen = int(row.get("intron_len", 0) or 0)
                except: ilen = 0
                arch = row.get("architecture", "?")
                akey, alab = arch_key(arch)
                meta = (f'<span class="arch-badge {akey}">{html.escape(alab)}</span>'
                        f'<span class="lbl">{html.escape(seqid)}:{start}-{end}({strand})</span>'
                        f'<span class="chip">len {ilen} bp</span>'
                        f'<span class="chip pur">telo {tfrac:.2f}</span>'
                        f'<span class="chip">bases {tbases}</span>'
                        f'<span class="chip">motif {html.escape(motif)}</span>')
                f.write(f'<div class="row"><div class="rmeta">{meta}</div>'
                        f'<pre class="seq">{colored}</pre></div>')

        f.write('</main></body></html>')
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
