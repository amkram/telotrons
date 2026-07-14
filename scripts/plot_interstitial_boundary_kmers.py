#!/usr/bin/env python3
"""Per-species figures: top 5'- and 3'-end k-mers of interstitial telomeric arrays.

Layout per species (one PNG):
    2 rows (5' end, 3' end) × 2 cols (6-mer, 12-mer)
    Each panel: top-N k-mers ranked by count

Admission criterion (Azzalin 2001 Chromosoma 110:75; Nergadze 2007;
Santagostino 2020 IJMS 21:2838; Giulotto-lab gold standard for ITS):
arrays must satisfy >=4 telomeric units AND <1 mismatch / unit
(identity to perfect tandem repeat >= 6/7). The raw
interstitial_arrays.tsv admits degenerate microsatellite noise via
find_arrays' fuzzy match; per MEMORY.md / its_detection_gold_standard,
73% of focal rows fail this criterion. Applied at top of main() to keep
the figure scripts in sync with interstitial_ortholog_textdump.py.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

from _common import slug as _slug, TEAL, RUST

TOP_N = 10

# ── ITS gold standard (Giulotto lab; mirrors interstitial_ortholog_textdump.py
#    MIN_ITS_UNITS / MIN_ITS_IDENT). Applied per row using the row's own
#    `motif` so it works for both TTAGGG (vertebrate / fungi) and TTTAGGG
#    (Apicomplexa / plants) panels.
MIN_ITS_UNITS = 4
MIN_ITS_IDENT = 6.0 / 7.0

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(s):
    return s.translate(_COMP)[::-1]


def _unit_identity(unit, motif):
    """Best identity of `unit` (one repeat-unit-length slice taken from the
    array boundary) to a perfect rotation of `motif` over fwd + revcomp.
    Returns 0.0 when inputs malformed so non-conforming rows are dropped."""
    if not unit or not motif:
        return 0.0
    u = unit.upper()
    m = motif.upper()
    if len(u) != len(m):
        return 0.0
    best = 0
    L = len(m)
    for base in (m, _rc(m)):
        for ph in range(L):
            rot = base[ph:] + base[:ph]
            hits = sum(1 for a, b in zip(u, rot) if a == b)
            if hits > best:
                best = hits
    return best / L


def _passes_its_gold_standard(row):
    """Giulotto criterion at the TSV-row level: array spans >= MIN_ITS_UNITS
    canonical units AND boundary units both match the motif at >= 6/7 identity.
    The boundary-unit check substitutes for whole-array tandem_identity (which
    requires the full array sequence, not exposed in interstitial_arrays.tsv);
    if either boundary unit is degenerate the row is treated as non-ITS noise."""
    motif = str(row.get("motif", "") or "")
    if not motif:
        return False
    try:
        alen = int(row.get("array_len", 0) or 0)
    except (TypeError, ValueError):
        return False
    if alen < MIN_ITS_UNITS * len(motif):
        return False
    fu = str(row.get("first_unit", "") or "")
    lu = str(row.get("last_unit", "") or "")
    if _unit_identity(fu, motif) < MIN_ITS_IDENT:
        return False
    if _unit_identity(lu, motif) < MIN_ITS_IDENT:
        return False
    return True


def plot_one(gid, organism, sub, outdir):
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), layout="constrained")
    n = len(sub)
    fig.suptitle(f"{organism}  ({gid})   {n} interstitial arrays  "
                 f"(non-genic, non-intronic, ≥5 kb from contig end)\n"
                 f"boundary k-mers extended 2 bp into flanking sequence",
                 fontsize=11, fontweight="bold")
    cells = [
        (0, 0, "first6",  "5' boundary  6 bp window  (array start −2 → +4)",  TEAL),
        (0, 1, "first12", "5' boundary  12 bp window  (array start −2 → +10)", TEAL),
        (1, 0, "last6",   "3' boundary  6 bp window  (array end −4 → +2)",   RUST),
        (1, 1, "last12",  "3' boundary  12 bp window  (array end −10 → +2)",  RUST),
    ]
    for r, c, col, title, color in cells:
        ax = axes[r, c]
        vc = sub[col].value_counts().head(TOP_N)
        if vc.empty:
            ax.text(0.5, 0.5, "no arrays", ha="center", va="center",
                    transform=ax.transAxes, color="#999")
            ax.set_xticks([]); ax.set_yticks([])
            for sp_ in ax.spines.values(): sp_.set_visible(False)
        else:
            labels = vc.index.tolist()[::-1]
            vals = vc.values.tolist()[::-1]
            ax.barh(labels, vals, color=color)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("array count")
            ax.tick_params(axis="y", labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)
    out = f"{outdir}/{_slug(gid)}_{_slug(organism)}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arrays", required=True, help="interstitial arrays TSV")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--filter-col", default=None)
    ap.add_argument("--filter-value", default=None)
    ap.add_argument("--min-records", type=int, default=5,
                    help="skip species with fewer arrays than this "
                         "(parity with plot_intron_boundary_kmers.py:--min-records=5 "
                         "and plot_interstitial_boundary_logos.py:--min-arrays=5)")
    ap.add_argument("--no-its-filter", action="store_true",
                    help="disable Giulotto >=4-unit / >=6/7 identity gold-standard "
                         "filter (default ON; matches interstitial_ortholog_textdump.py)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.arrays, sep="\t", low_memory=False)
    if args.filter_col and args.filter_value is not None:
        before = len(df)
        vals = [v for v in str(args.filter_value).split(",") if v]
        df = df[df[args.filter_col].astype(str).isin(vals)]
        print(f"filter {args.filter_col} in {vals}: {before} → {len(df)}",
              flush=True)
    # ITS gold-standard filter (Giulotto >=4 units, >=6/7 identity).
    if not args.no_its_filter:
        before = len(df)
        keep = df.apply(_passes_its_gold_standard, axis=1)
        df = df[keep]
        print(f"ITS gold-standard filter (>=4 units, >=6/7 identity): "
              f"{before} → {len(df)}", flush=True)
    if df.empty:
        print("no interstitial arrays; nothing to plot")
        return
    n_written = 0
    for gid, sub in df.groupby("genome_id"):
        organism = sub.organism.iloc[0]
        if len(sub) < args.min_records:
            print(f"  {gid}: {len(sub)} < --min-records={args.min_records}; skip")
            continue
        plot_one(gid, organism, sub, args.outdir)
        print(f"  {gid}: {len(sub)}")
        n_written += 1
    print(f"wrote {n_written} figures into {args.outdir}/")


if __name__ == "__main__":
    main()
