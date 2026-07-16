#!/usr/bin/env python3
"""Render an HTML preview of the current telotron candidate set for inspection.

Per species, sorted by (telomeric_frac desc, telomeric_bases desc, intron_len desc)
so the longest and most-telomeric introns land at the top of each species panel.
50 bp flanks either side of the intron.

Color scheme (motif rotations + reverse complements, exact + 1-mismatch):
    forest green      G-rich exact motif rotation
    pale green        G-rich 1-mismatch rotation
    crimson           C-rich exact motif rotation (reverse-complement)
    pale red          C-rich 1-mismatch rotation
    mid-grey          flank (non-motif)
    dark grey         intron (non-motif)
Intron boundaries marked with a thin `⎢` divider.
"""
import argparse
import csv
import gzip
import html
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import find_genome_fasta  # noqa: E402

FLANK = 50
RC_TAB = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(RC_TAB)[::-1]


def rotations(m: str) -> set:
    return {m[i:] + m[:i] for i in range(len(m))}


def load_fasta(path: str) -> dict:
    """Load a FASTA (plain or .gz) into {seqid: uppercase_seq}."""
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
    """Return (fwd_exact, fwd_1mm, rev_exact, rev_1mm) sets of kmers of length len(motif)."""
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

    fwd_1mm = mismatches(fwd) - fwd - rev
    rev_1mm = mismatches(rev) - rev - fwd
    return fwd, fwd_1mm, rev, rev_1mm


CLS_RANK = {"ge": 5, "ce": 4, "gm": 3, "cm": 2}


def classify_positions(seq: str, motif: str):
    """Per-position class: 'ge' (G exact), 'ce' (C exact), 'gm' (G 1mm), 'cm' (C 1mm), or None."""
    n = len(seq)
    L = len(motif)
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
    """HTML colored spans with intron-boundary markers."""
    n = len(seq)
    intron_start = up_len
    intron_end = n - dn_len
    cls = classify_positions(seq, motif)
    # fill non-motif positions with 'fl' or 'in' by region
    eff = []
    for i in range(n):
        c = cls[i]
        if c is None:
            c = "fl" if (i < intron_start or i >= intron_end) else "in"
        eff.append(c)
    parts = []
    buf, cur = [], None
    for i in range(n):
        # Insert boundary marker before intron_start and before intron_end (so it sits between adjacent bases)
        if i == intron_start:
            if buf:
                parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
                buf, cur = [], None
            parts.append('<span class="sep">⎢</span>')
        if i == intron_end and intron_end > intron_start:
            if buf:
                parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
                buf, cur = [], None
            parts.append('<span class="sep">⎢</span>')
        c = eff[i]
        if c != cur:
            if buf:
                parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
                buf = []
            cur = c
        buf.append(html.escape(seq[i]))
    if buf:
        parts.append(f'<span class="{cur}">{"".join(buf)}</span>')
    return "".join(parts)


CSS = """
:root {
  --bg: #fbfaf7;
  --fg: #16161a;
  --muted: #6a6a72;
  --card: #ffffff;
  --line: #e6e4dd;
  --accent: #1b5f3a;
  --g-exact: #1b7a3a;
  --g-mm:    #7fb992;
  --c-exact: #b2182b;
  --c-mm:    #d68b9a;
  --flank:   #999;
  --intron:  #2b2b2f;
  --sep:     #a89f7c;
  --chip-bg: #efece5;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #17181b; --fg: #ececee; --muted: #8a8a91;
    --card: #1e2024; --line: #2a2c31; --accent: #7ec89d;
    --g-exact: #5ec97e; --g-mm: #a9d5b9;
    --c-exact: #ff7683; --c-mm: #f0b6bd;
    --flank: #7d7d85; --intron: #c9c9cd; --sep: #c2b477;
    --chip-bg: #262930;
  }
}
:root[data-theme="light"] {
  --bg: #fbfaf7; --fg: #16161a; --muted: #6a6a72;
  --card: #ffffff; --line: #e6e4dd; --accent: #1b5f3a;
  --g-exact: #1b7a3a; --g-mm: #7fb992;
  --c-exact: #b2182b; --c-mm: #d68b9a;
  --flank: #999; --intron: #2b2b2f; --sep: #a89f7c;
  --chip-bg: #efece5;
}
:root[data-theme="dark"] {
  --bg: #17181b; --fg: #ececee; --muted: #8a8a91;
  --card: #1e2024; --line: #2a2c31; --accent: #7ec89d;
  --g-exact: #5ec97e; --g-mm: #a9d5b9;
  --c-exact: #ff7683; --c-mm: #f0b6bd;
  --flank: #7d7d85; --intron: #c9c9cd; --sep: #c2b477;
  --chip-bg: #262930;
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: var(--fg);
  font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
  font-size: 14px; line-height: 1.5;
}
header {
  padding: 24px 40px 16px;
  border-bottom: 1px solid var(--line);
  background: var(--card);
  position: sticky; top: 0; z-index: 10;
  backdrop-filter: saturate(180%) blur(6px);
}
h1 { font-size: 20px; font-weight: 600; margin: 0 0 4px; letter-spacing: -0.01em; }
h2 {
  font-size: 15px; font-weight: 600; margin: 28px 0 4px;
  display: flex; align-items: baseline; gap: 12px;
  padding-top: 8px;
}
.org { font-weight: 400; font-style: italic; color: var(--muted); font-size: 13px; }
p.meta { color: var(--muted); font-size: 12.5px; margin: 4px 0 6px; }
.legend {
  display: flex; gap: 16px; align-items: center; margin-top: 8px;
  font-size: 11.5px; color: var(--muted);
  font-family: ui-monospace, monospace;
}
.legend .sw { display: inline-block; width: 12px; height: 12px;
              vertical-align: -2px; margin-right: 4px; border-radius: 2px; }
.jump {
  display: flex; flex-wrap: wrap; gap: 4px 6px; margin-top: 12px;
  max-height: 96px; overflow-y: auto;
}
.jump a {
  color: var(--fg); text-decoration: none;
  padding: 3px 8px; border: 1px solid var(--line); border-radius: 4px;
  font-size: 11px; display: inline-flex; align-items: center; gap: 6px;
  background: var(--bg);
}
.jump a:hover { border-color: var(--accent); }
.jump b { font-weight: 500; font-family: ui-monospace, monospace; }
.jump .ct { color: var(--muted); font-variant-numeric: tabular-nums; }
.jump .src {
  color: var(--muted); font-size: 10px;
  text-transform: uppercase; letter-spacing: 0.05em;
}
section { padding: 4px 40px 20px; }
.row {
  padding: 8px 0;
  border-bottom: 1px solid var(--line);
}
.rmeta {
  font-size: 11.5px; display: flex; flex-wrap: wrap;
  gap: 4px 10px; align-items: baseline;
}
.lbl { font-family: ui-monospace, monospace; color: var(--muted); }
.chip {
  background: var(--chip-bg); padding: 1px 6px; border-radius: 3px;
  font-size: 10.5px; color: var(--muted);
  font-variant-numeric: tabular-nums;
}
.chip.arch { color: var(--accent); }
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
"""


def build_html(rows_by_species, canonical, species_order, out_path, refseq_dir, tara_dir, max_per_species):
    fasta_cache = {}

    def get_seq(gid, seqid):
        if gid not in fasta_cache:
            p = find_genome_fasta(gid, refseq_dir, tara_dir, required=False)
            fasta_cache[gid] = load_fasta(p) if p else {}
        return fasta_cache[gid].get(seqid)

    total = sum(len(v) for v in rows_by_species.values())
    with open(out_path, "w") as f:
        f.write('<meta charset="utf-8">'
                '<title>Telotron candidate preview</title>'
                f'<style>{CSS}</style></head><body>')
        f.write('<header>')
        f.write('<h1>Telotron candidate set — visual inspection</h1>')
        f.write(f'<p class="meta">{total} loci across {len(rows_by_species)} species '
                'sorted per-species by telomeric_frac desc, then telomeric_bases desc, then intron_len desc. '
                '±50 bp flanks. Intron boundaries marked with <span style="color:var(--sep);font-family:ui-monospace,monospace">⎢</span>.</p>')
        f.write('<div class="legend">'
                '<span><span class="sw" style="background:var(--g-exact)"></span>G-rich exact</span>'
                '<span><span class="sw" style="background:var(--g-mm)"></span>G-rich 1-mm</span>'
                '<span><span class="sw" style="background:var(--c-exact)"></span>C-rich exact</span>'
                '<span><span class="sw" style="background:var(--c-mm)"></span>C-rich 1-mm</span>'
                '<span><span class="sw" style="background:var(--flank)"></span>flank</span>'
                '<span><span class="sw" style="background:var(--intron)"></span>intron (non-motif)</span>'
                '</div>')
        f.write('<nav class="jump">')
        for gid in species_order:
            org = rows_by_species[gid][0].get("organism", "")
            n = len(rows_by_species[gid])
            src = rows_by_species[gid][0].get("source", "?")
            slug = gid.replace(".", "_").replace("/", "_")
            f.write(f'<a href="#s-{slug}"><b>{html.escape(gid)}</b>'
                    f'<span class="ct">{n}</span>'
                    f'<span class="src">{src}</span></a>')
        f.write('</nav></header>')

        for gid in species_order:
            slug = gid.replace(".", "_").replace("/", "_")
            org = rows_by_species[gid][0].get("organism", "")
            src = rows_by_species[gid][0].get("source", "?")
            motif_c = canonical.get(gid, "") or ""
            all_rows = rows_by_species[gid]
            n = len(all_rows)
            display = all_rows[:max_per_species]
            capped = n > max_per_species
            f.write(f'<section id="s-{slug}">'
                    f'<h2>{html.escape(gid)} <span class="org">{html.escape(org)}</span></h2>')
            note = f' (top {max_per_species} shown)' if capped else ''
            canon_html = html.escape(motif_c) if motif_c else '<em>(none — per-locus motif from row)</em>'
            f.write(f'<p class="meta">{src} · {n} loci{note} · canonical motif: <code>{canon_html}</code></p>')

            for r in display:
                try:
                    start, end = int(r["start"]), int(r["end"])
                except (KeyError, ValueError):
                    continue
                seqid = r["seqid"]
                full = get_seq(gid, seqid)
                if not full:
                    f.write(f'<div class="row missing">{html.escape(seqid)}:{start}-{end} — seq not found</div>')
                    continue
                # GFF is 1-based inclusive; python slicing is 0-based half-open.
                # Intron occupies python[start-1 : end]. Boundary markers are the
                # whole point of this preview, so being one base off silently
                # puts the GT donor into the left-flank block.
                pystart, pyend = start - 1, end
                up_s = max(0, pystart - FLANK)
                dn_e = min(len(full), pyend + FLANK)
                seg = full[up_s:dn_e]
                up_len = pystart - up_s
                dn_len = dn_e - pyend
                strand = r.get("strand", "+")
                if strand == "-":
                    seg = rc(seg)
                    up_len, dn_len = dn_len, up_len
                motif = (r.get("motif") or motif_c or "TTAGGG").strip()
                colored = render_seq(seg, motif, up_len, dn_len)
                try:
                    tfrac = float(r["telomeric_frac"])
                except (KeyError, ValueError):
                    tfrac = 0.0
                try:
                    tbases = int(float(r.get("telomeric_bases", 0)))
                except (ValueError, TypeError):
                    tbases = 0
                try:
                    ilen = int(r.get("intron_len", 0) or 0)
                except (ValueError, TypeError):
                    ilen = 0
                meta = (f'<span class="lbl">{html.escape(seqid)}:{start}-{end}({strand})</span>'
                        f'<span class="chip">len {ilen} bp</span>'
                        f'<span class="chip pur">telo {tfrac:.2f}</span>'
                        f'<span class="chip">bases {tbases}</span>'
                        f'<span class="chip">motif {html.escape(motif)}</span>'
                        f'<span class="chip arch">{html.escape(r.get("admission_pathway","?"))}</span>')
                f.write(f'<div class="row"><div class="rmeta">{meta}</div>'
                        f'<pre class="seq">{colored}</pre></div>')
            f.write('</section>')
        f.write('</body></html>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", default="work/results/final_telotron_set.tsv")
    ap.add_argument("--refseq-dir", default="data/raw/refseq")
    ap.add_argument("--tara-dir", default="data/raw/tara")
    ap.add_argument("--canonical", default="work/manifests/canonical_motifs.tsv")
    ap.add_argument("--out", default="work/results/candidate_preview.html")
    ap.add_argument("--max-per-species", type=int, default=200)
    args = ap.parse_args()

    canonical = {}
    if os.path.exists(args.canonical):
        for r in csv.DictReader(open(args.canonical), delimiter="\t"):
            if r.get("motif"):
                canonical[r["genome_id"]] = r["motif"]

    rows_by = defaultdict(list)
    with open(args.final) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows_by[r["genome_id"]].append(r)

    species_order = sorted(rows_by, key=lambda g: -len(rows_by[g]))
    for gid in species_order:
        rows_by[gid].sort(key=lambda r: (
            -float(r.get("telomeric_frac", 0) or 0),
            -int(float(r.get("telomeric_bases", 0) or 0)),
            -int(r.get("intron_len", 0) or 0),
        ))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    build_html(rows_by, canonical, species_order, args.out,
               args.refseq_dir, args.tara_dir, args.max_per_species)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
