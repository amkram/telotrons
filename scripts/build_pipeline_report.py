#!/usr/bin/env python3
"""Single self-contained HTML report of the telotrons pipeline.

Sections (each <details>, collapsed by default except Overview):
    1. Overview                   counts, species breakdown, pipeline metadata
    2. Final telotron set         sortable table
    3. Per-species summary        sortable table
    4. Architecture               per-architecture counts + figure
    5. Boundary k-mer enrichment  table + figures
    6. Distance to contig end     table + figure
    7. Pipeline stages            stage figures
    8. Composite per-species      composite figures (kmers/logos × all/5p3p)
    9. MSAs                       per-(species, arch) combined alignment colored
   10. Linker BLAST hits          own-genome + all-genomes summaries
   11. Ortholog locus_text        per-locus collapsible blocks (sample)
   12. Telogator2 (long read)     per-species telomere lengths if present

Embeds PNGs as base64 data URIs so the HTML is fully self-contained and
shareable. Large tables (>500 rows) are truncated with a "show all" link
or skipped depending on `--max-table-rows`.
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import os
import re
import sys
from pathlib import Path

NUC_COLOR = {
    "A": "#7ec97e",  # green
    "C": "#7eb6ff",  # blue
    "G": "#ffc97e",  # orange
    "T": "#ff8a8a",  # red
    "U": "#ff8a8a",
    "N": "#cccccc",
    "-": "#eeeeee",
    " ": "transparent",
}


CSS = """
* { box-sizing: border-box; }
body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    margin: 0;
    padding: 0 24px 80px;
    color: #1f2933;
    background: #fafafa;
    line-height: 1.5;
}
header {
    background: #20808d;
    color: white;
    padding: 18px 24px;
    margin: 0 -24px 24px;
    border-bottom: 2px solid #186370;
}
header h1 { margin: 0; font-size: 22px; font-weight: 600; }
header .meta { font-size: 12px; opacity: 0.85; margin-top: 4px; }
nav {
    position: sticky; top: 0;
    background: white;
    margin: 0 -24px;
    padding: 6px 24px;
    border-bottom: 1px solid #e0e6ec;
    z-index: 100;
    font-size: 13px;
}
nav a { color: #20808d; text-decoration: none; margin-right: 14px; }
nav a:hover { text-decoration: underline; }
details {
    margin: 16px 0;
    background: white;
    border: 1px solid #e0e6ec;
    border-radius: 6px;
}
details > summary {
    cursor: pointer;
    padding: 12px 16px;
    font-weight: 600;
    font-size: 16px;
    list-style: none;
    user-select: none;
    background: #f4f7fa;
    border-radius: 6px 6px 0 0;
}
details[open] > summary { border-bottom: 1px solid #e0e6ec; }
details > summary::before {
    content: "▶";
    display: inline-block;
    transition: transform 0.15s;
    margin-right: 8px;
    color: #20808d;
    font-size: 11px;
}
details[open] > summary::before { transform: rotate(90deg); }
.body { padding: 14px 18px; }
table {
    border-collapse: collapse;
    width: 100%;
    font-size: 13px;
    margin-top: 8px;
}
th, td {
    border-bottom: 1px solid #e8ecef;
    padding: 6px 10px;
    text-align: left;
}
th {
    background: #f4f7fa;
    cursor: pointer;
    user-select: none;
    position: sticky;
    top: 0;
}
th.sorted-asc::after { content: " ↑"; }
th.sorted-desc::after { content: " ↓"; }
tr:hover td { background: #fafbfc; }
.table-wrap {
    max-height: 600px;
    overflow: auto;
    border: 1px solid #e0e6ec;
    border-radius: 4px;
}
img.figure {
    max-width: 100%;
    height: auto;
    border: 1px solid #e0e6ec;
    border-radius: 4px;
    margin-top: 8px;
}
.figure-row {
    display: flex;
    flex-wrap: wrap;
    gap: 16px;
    margin-top: 8px;
}
.figure-row > div { flex: 1 1 360px; }
.figure-row figcaption {
    font-size: 12px;
    color: #5d6975;
    margin-bottom: 4px;
}
pre.aln {
    background: #fafafa;
    border: 1px solid #e0e6ec;
    padding: 8px;
    font-family: "SF Mono", "Cascadia Mono", Menlo, Consolas, monospace;
    font-size: 11px;
    line-height: 1.25;
    overflow-x: auto;
    white-space: pre;
}
pre.aln span { padding: 0 0.5px; }
.locus-text {
    background: #fafafa;
    border: 1px solid #e0e6ec;
    padding: 10px;
    font-family: "SF Mono", "Cascadia Mono", Menlo, Consolas, monospace;
    font-size: 11px;
    white-space: pre-wrap;
    overflow-x: auto;
    max-height: 480px;
}
.stat-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 12px;
    margin: 12px 0;
}
.stat {
    background: #f4f7fa;
    padding: 12px 14px;
    border-radius: 4px;
    border-left: 3px solid #20808d;
}
.stat .label { font-size: 11px; color: #5d6975; text-transform: uppercase; letter-spacing: 0.04em; }
.stat .value { font-size: 22px; font-weight: 600; margin-top: 2px; }
.muted { color: #5d6975; font-size: 12px; }
.note {
    background: #fef9e7;
    border-left: 3px solid #f0b429;
    padding: 8px 12px;
    margin: 8px 0;
    font-size: 12px;
}
nav .arch-key { font-family: monospace; }
.legend {
    display: inline-flex; gap: 10px; align-items: center;
    font-size: 11px; color: #5d6975; margin: 4px 0 8px;
}
.legend span.swatch {
    display: inline-block; width: 14px; height: 14px;
    border: 1px solid #e0e6ec; vertical-align: middle;
}
"""

JS = r"""
// Sortable table: click <th> to sort.
document.querySelectorAll("table.sortable").forEach((table) => {
  const headers = table.querySelectorAll("thead th");
  headers.forEach((th, idx) => {
    th.addEventListener("click", () => {
      const tbody = table.querySelector("tbody");
      const rows = Array.from(tbody.querySelectorAll("tr"));
      const asc = !th.classList.contains("sorted-asc");
      headers.forEach((h) => h.classList.remove("sorted-asc", "sorted-desc"));
      th.classList.add(asc ? "sorted-asc" : "sorted-desc");
      rows.sort((a, b) => {
        const av = a.children[idx].textContent.trim();
        const bv = b.children[idx].textContent.trim();
        const an = parseFloat(av);
        const bn = parseFloat(bv);
        if (!isNaN(an) && !isNaN(bn)) return asc ? an - bn : bn - an;
        return asc ? av.localeCompare(bv) : bv.localeCompare(av);
      });
      rows.forEach((r) => tbody.appendChild(r));
    });
  });
});
"""


# ──────────────────────────────────────────────────────────────────────────
# Generic helpers
# ──────────────────────────────────────────────────────────────────────────

def _read_tsv(path: Path, max_rows: int | None = None):
    if not path.exists() or path.stat().st_size == 0:
        return None, []
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return None, []
        rows = []
        for i, row in enumerate(reader):
            if max_rows is not None and i >= max_rows:
                break
            rows.append(row)
    return header, rows


def _tsv_total_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as fh:
        return max(0, sum(1 for _ in fh) - 1)


def _img_b64(path: Path) -> str | None:
    if not path.exists():
        return None
    mime = "image/png" if path.suffix.lower() == ".png" else \
           "image/jpeg" if path.suffix.lower() in (".jpg", ".jpeg") else \
           "image/svg+xml" if path.suffix.lower() == ".svg" else \
           "image/png"
    data = base64.b64encode(path.read_bytes()).decode()
    return f"data:{mime};base64,{data}"


# Globals patched by main() so renderers can read them without passing args around.
_OPTS = {"embed_images": False, "out_dir": Path(".")}


def _img_tag(path: Path, alt: str = "", caption: str | None = None) -> str:
    if not path.exists():
        return f'<div class="muted">[missing: {html.escape(str(path))}]</div>'
    if _OPTS["embed_images"]:
        src = _img_b64(path) or ""
    else:
        try:
            src = os.path.relpath(path, _OPTS["out_dir"])
        except ValueError:
            src = str(path)
    cap_html = f'<figcaption>{html.escape(caption)}</figcaption>' if caption else ""
    return (
        f'<div><figure style="margin:0">{cap_html}'
        f'<img class="figure" loading="lazy" alt="{html.escape(alt)}" src="{html.escape(src)}"/></figure></div>'
    )


def _table_html(header, rows, table_id: str | None = None, sortable: bool = True) -> str:
    if not header:
        return '<div class="muted">[no data]</div>'
    cls = "sortable" if sortable else ""
    id_attr = f' id="{table_id}"' if table_id else ""
    thead = "".join(f"<th>{html.escape(c)}</th>" for c in header)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in row) + "</tr>"
        for row in rows
    )
    return (
        f'<div class="table-wrap">'
        f'<table class="{cls}"{id_attr}>'
        f'<thead><tr>{thead}</tr></thead>'
        f'<tbody>{body}</tbody>'
        f'</table>'
        f'</div>'
    )


def _section(title: str, body: str, open: bool = False, anchor: str | None = None) -> str:
    anchor_attr = f' id="{anchor}"' if anchor else ""
    open_attr = " open" if open else ""
    return f'<details{open_attr}{anchor_attr}><summary>{html.escape(title)}</summary><div class="body">{body}</div></details>'


def _color_seq(seq: str) -> str:
    """Render a nucleotide sequence with per-base coloring."""
    out = []
    for c in seq:
        up = c.upper()
        color = NUC_COLOR.get(up, "transparent")
        if color == "transparent":
            out.append(html.escape(c))
        else:
            out.append(f'<span style="background:{color}">{html.escape(c)}</span>')
    return "".join(out)


# ──────────────────────────────────────────────────────────────────────────
# Section renderers
# ──────────────────────────────────────────────────────────────────────────

def section_overview(results: Path) -> str:
    species_tsv = results / "final_species_summary.tsv"
    final_tsv = results / "final_telotron_set.tsv"
    arch_tsv = results / "architecture_summary.tsv"

    species_n = _tsv_total_rows(species_tsv)
    final_n = _tsv_total_rows(final_tsv)
    pos_species_n = 0
    if species_tsv.exists():
        header, rows = _read_tsv(species_tsv)
        if header and "telotrons" in header:
            tcol = header.index("telotrons")
            pos_species_n = sum(1 for r in rows if r[tcol] not in ("", "0", "0.0"))

    arch_counts = {}
    if arch_tsv.exists():
        header, rows = _read_tsv(arch_tsv)
        if header and "architecture" in header and "n_loci" in header:
            aidx = header.index("architecture")
            nidx = header.index("n_loci")
            for r in rows:
                arch_counts[r[aidx]] = arch_counts.get(r[aidx], 0) + int(r[nidx] or 0)

    stat_html = (
        f'<div class="stat-grid">'
        f'<div class="stat"><div class="label">Species in manifest</div><div class="value">{species_n}</div></div>'
        f'<div class="stat"><div class="label">Positive species</div><div class="value">{pos_species_n}</div></div>'
        f'<div class="stat"><div class="label">Final telotron loci</div><div class="value">{final_n}</div></div>'
        f'<div class="stat"><div class="label">Pipeline outputs</div><div class="value">{sum(1 for _ in results.rglob("*.tsv"))}</div><div class="muted">TSVs</div></div>'
        f'</div>'
    )

    arch_html = ""
    if arch_counts:
        rows = "".join(
            f"<tr><td>{html.escape(a)}</td><td>{n}</td></tr>"
            for a, n in sorted(arch_counts.items(), key=lambda x: -x[1])
        )
        arch_html = (
            f'<h3>Architecture breakdown</h3>'
            f'<div class="table-wrap"><table class="sortable">'
            f'<thead><tr><th>architecture</th><th>n_loci</th></tr></thead>'
            f'<tbody>{rows}</tbody></table></div>'
        )

    return stat_html + arch_html


def section_tsv(label: str, results: Path, name: str, max_rows: int = 500) -> str:
    path = results / name
    if not path.exists():
        return f'<div class="muted">[missing: {name}]</div>'
    total = _tsv_total_rows(path)
    header, rows = _read_tsv(path, max_rows=max_rows)
    note = ""
    if total > max_rows:
        note = f'<div class="note">{label}: showing first {max_rows} of {total} rows. Full TSV at {html.escape(name)}.</div>'
    return note + _table_html(header, rows)


def section_figures(results: Path, glob: str, label: str = "", max_figs: int = 60) -> str:
    figs = sorted((results / "figures").glob(glob))
    if not figs:
        return f'<div class="muted">[no figures matched {glob}]</div>'
    truncated = len(figs) > max_figs
    figs_to_show = figs[:max_figs]
    items = [_img_tag(f, alt=f.stem, caption=f.name) for f in figs_to_show]
    note = ""
    if truncated:
        note = (f'<div class="note">Showing first {max_figs} of {len(figs)} figures. '
                f'Full glob: {html.escape(glob)}</div>')
    return note + f'<div class="figure-row">{"".join(items)}</div>'


def _color_combined_line(line: str) -> str:
    """A combined.aln.* data line is `LXXXX  seq1 seq2 seq3 ...`.
    Color the sequence portions, leave the locus ID + spaces uncolored.
    """
    if not line:
        return ""
    parts = line.split("  ", 1)
    if len(parts) == 1:
        return _color_seq(line)
    locus_id, rest = parts
    return html.escape(locus_id) + "  " + _color_seq(rest)


def section_msa(results: Path, max_loci: int = 6) -> str:
    """Render combined.aln.{txt,fa} with nucleotide coloring.

    The combined.aln.* file is a column-aligned plain-text view with
    region boundaries shown as single-space gaps. We color the data lines and
    keep the header (#-prefixed lines) verbatim.
    """
    msa_root = results / "msa_regions"
    if not msa_root.exists():
        return '<div class="muted">[no msa_regions/ found]</div>'
    out = []
    legend = (
        '<div class="legend">'
        f'<span class="swatch" style="background:{NUC_COLOR["A"]}"></span>A '
        f'<span class="swatch" style="background:{NUC_COLOR["C"]}"></span>C '
        f'<span class="swatch" style="background:{NUC_COLOR["G"]}"></span>G '
        f'<span class="swatch" style="background:{NUC_COLOR["T"]}"></span>T '
        f'<span class="swatch" style="background:{NUC_COLOR["-"]}"></span>gap '
        f'</div>'
    )
    out.append(legend)
    for species_dir in sorted(msa_root.iterdir()):
        if not species_dir.is_dir():
            continue
        species_blocks = []
        for arch_dir in sorted(species_dir.iterdir()):
            if not arch_dir.is_dir():
                continue
            combined = next(
                (arch_dir / n for n in ("combined.aln.txt", "combined.aln.fa")
                 if (arch_dir / n).exists()),
                None,
            )
            if combined is None:
                continue
            txt = combined.read_text()
            lines = txt.splitlines()
            header_lines = [l for l in lines if l.startswith("#")]
            data_lines = [l for l in lines if not l.startswith("#") and l.strip()]
            data_lines = data_lines[:max_loci]
            if not data_lines:
                continue
            colored = "\n".join(_color_combined_line(l) for l in data_lines)
            header_block = "\n".join(html.escape(h) for h in header_lines)
            species_blocks.append(
                f'<details><summary>{html.escape(arch_dir.name)} '
                f'({len(data_lines)} loci shown)</summary>'
                f'<pre class="aln">{header_block}\n{colored}</pre></details>'
            )
        if species_blocks:
            out.append(
                f'<details><summary>{html.escape(species_dir.name)} '
                f'({len(species_blocks)} architectures)</summary>'
                f'<div class="body">{"".join(species_blocks)}</div></details>'
            )
    if len(out) == 1:
        return out[0] + '<div class="muted">[no MSA outputs]</div>'
    return "".join(out)


def section_locus_text(results: Path, max_files: int = 20) -> str:
    """Sample the v2 ortholog locus_text output."""
    root = results / "telotron_orthologs_v2" / "locus_text"
    if not root.exists():
        return '<div class="muted">[no telotron_orthologs_v2/locus_text/ found]</div>'
    files = sorted(root.rglob("*.txt"))[:max_files]
    if not files:
        return '<div class="muted">[no .txt files under locus_text]</div>'
    out = [f'<div class="muted">Showing first {len(files)} of '
           f'{sum(1 for _ in root.rglob("*.txt"))} locus_text files.</div>']
    for f in files:
        rel = f.relative_to(root)
        body = f.read_text()
        out.append(
            f'<details><summary>{html.escape(str(rel))}</summary>'
            f'<pre class="locus-text">{html.escape(body)}</pre></details>'
        )
    return "".join(out)


def section_blast_summary(results: Path) -> str:
    own = results / "linker_blast_hits_own_genome.tsv"
    allgen = results / "linker_blast_hits_all_genomes.tsv"
    parts = []
    for label, path in (("Own-genome hits", own), ("Whole-DB hits", allgen)):
        if not path.exists():
            continue
        total = _tsv_total_rows(path)
        header, rows = _read_tsv(path, max_rows=300)
        parts.append(
            f'<h3>{html.escape(label)}</h3>'
            f'<div class="muted">{total:,} rows total. First 300 shown.</div>'
            + _table_html(header, rows)
        )
    return "".join(parts) or '<div class="muted">[no linker BLAST results]</div>'


def section_telogator2(results: Path) -> str:
    tg2 = results / "telogator2"
    if not tg2.exists():
        return '<div class="muted">[no telogator2/ output]</div>'
    parts = []
    for run_dir in sorted(tg2.iterdir()):
        if not run_dir.is_dir():
            continue
        tlens = run_dir / "tlens_by_allele.tsv"
        if not tlens.exists():
            continue
        header, rows = _read_tsv(tlens, max_rows=200)
        total = _tsv_total_rows(tlens)
        parts.append(
            f'<h3>{html.escape(run_dir.name)}</h3>'
            f'<div class="muted">{total:,} alleles total. First 200 shown.</div>'
            + _table_html(header, rows)
        )
    return "".join(parts) or '<div class="muted">[no telomere length tables]</div>'


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="work/results", help="pipeline results dir")
    ap.add_argument("--out", default="work/results/pipeline_report.html")
    ap.add_argument("--max-table-rows", type=int, default=500)
    ap.add_argument("--max-msa-loci", type=int, default=6)
    ap.add_argument("--max-figs-per-section", type=int, default=60)
    ap.add_argument("--embed-images", action="store_true",
                    help="base64-embed every figure (self-contained HTML, large file). "
                         "Default: link via relative paths.")
    args = ap.parse_args()

    results = Path(args.results)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    _OPTS["embed_images"] = args.embed_images
    _OPTS["out_dir"] = out.parent

    sections = []

    # Nav anchors
    nav_html = (
        '<nav>'
        '<a href="#overview">Overview</a>'
        '<a href="#final-set">Final set</a>'
        '<a href="#species">Per-species</a>'
        '<a href="#architecture">Architecture</a>'
        '<a href="#kmers">Boundary k-mers</a>'
        '<a href="#distance">Distance</a>'
        '<a href="#stages">Pipeline stages</a>'
        '<a href="#composites">Composites</a>'
        '<a href="#msa">MSAs</a>'
        '<a href="#blast">Linker BLAST</a>'
        '<a href="#orthologs">Orthologs</a>'
        '<a href="#telogator">Telogator2</a>'
        '</nav>'
    )

    sections.append(_section(
        "1. Overview", section_overview(results), open=True, anchor="overview",
    ))
    sections.append(_section(
        "2. Final telotron set",
        section_tsv("Final telotron set", results, "final_telotron_set.tsv",
                    max_rows=args.max_table_rows),
        anchor="final-set",
    ))
    sections.append(_section(
        "3. Per-species summary",
        section_tsv("Per-species summary", results, "final_species_summary.tsv",
                    max_rows=args.max_table_rows),
        anchor="species",
    ))
    sections.append(_section(
        "4. Architecture",
        section_tsv("Architecture summary", results, "architecture_summary.tsv",
                    max_rows=args.max_table_rows)
        + section_figures(results, "boundary_kmers_by_arch/*.png", max_figs=args.max_figs_per_section),
        anchor="architecture",
    ))
    sections.append(_section(
        "5. Boundary k-mer enrichment",
        section_tsv("Boundary k-mer enrichment", results, "boundary_kmer_enrichment.tsv",
                    max_rows=args.max_table_rows)
        + section_figures(results, "telotron_boundary_kmers/*.png", label="telotron_boundary_kmers", max_figs=args.max_figs_per_section)
        + section_figures(results, "non_telotron_boundary_kmers/*.png", label="non_telotron_boundary_kmers", max_figs=args.max_figs_per_section),
        anchor="kmers",
    ))
    sections.append(_section(
        "6. Distance to contig end",
        section_tsv("Distance to end", results, "distance_to_end.tsv",
                    max_rows=args.max_table_rows),
        anchor="distance",
    ))
    sections.append(_section(
        "7. Pipeline stages",
        section_figures(results, "pipeline_stages/*.png", max_figs=args.max_figs_per_section),
        anchor="stages",
    ))
    sections.append(_section(
        "8. Composite per-species",
        section_figures(results, "composite_boundary_kmers*/*.png", max_figs=args.max_figs_per_section)
        + section_figures(results, "composite_boundary_logos*/*.png", max_figs=args.max_figs_per_section),
        anchor="composites",
    ))
    sections.append(_section(
        "9. MSAs (colored alignments)",
        section_msa(results, max_loci=args.max_msa_loci),
        anchor="msa",
    ))
    sections.append(_section(
        "10. Linker BLAST hits",
        section_blast_summary(results),
        anchor="blast",
    ))
    sections.append(_section(
        "11. Ortholog locus_text (sample)",
        section_locus_text(results),
        anchor="orthologs",
    ))
    sections.append(_section(
        "12. Telogator2 (long-read telomere lengths)",
        section_telogator2(results),
        anchor="telogator",
    ))

    header_html = (
        '<header>'
        '<h1>Telotrons — pipeline report</h1>'
        f'<div class="meta">generated from {html.escape(str(results))}</div>'
        '</header>'
    )

    html_doc = (
        '<!DOCTYPE html><html lang="en"><head>'
        '<meta charset="utf-8"/>'
        '<meta name="viewport" content="width=device-width, initial-scale=1"/>'
        '<title>Telotrons pipeline report</title>'
        f'<style>{CSS}</style>'
        '</head><body>'
        + header_html
        + nav_html
        + "".join(sections)
        + f'<script>{JS}</script>'
        '</body></html>'
    )

    out.write_text(html_doc)
    size_mb = out.stat().st_size / 1024 / 1024
    print(f"wrote {out} ({size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
