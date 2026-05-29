#!/usr/bin/env python3
"""Per-species figures: top boundary 6-mers / 12-mers at intron ends.

Generic over any intron-style TSV that has `genome_id`, `organism`, `first40`,
`last40` columns — works on the telotron set, the non-telotron control set,
or any subset. Layout mirrors plot_interstitial_boundary_kmers.py:

    2 rows (donor, acceptor) × 2 cols (6-mer, 12-mer)
    Each panel: top-N k-mers ranked by count

Boundary windows are sliced straight from `first40` (donor: prefix bases,
spans GT and downstream) and `last40` (acceptor: suffix bases, ends in AG).
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

from _common import slug as _slug, TEAL, RUST

TOP_N = 10


def plot_one(label, sub, outdir, title_suffix="", out_name=None, title=None):
    n = len(sub)
    if n == 0:
        return False
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), layout="constrained")
    fig.suptitle(f"{title or label}   n={n}{title_suffix}",
                 fontsize=11, fontweight="bold")
    donor6   = sub.first40.astype(str).str.slice(0, 6)
    donor12  = sub.first40.astype(str).str.slice(0, 12)
    accept6  = sub.last40.astype(str).str.slice(-6)
    accept12 = sub.last40.astype(str).str.slice(-12)
    cells = [
        (0, 0, donor6,   "5' donor  6 bp  (intron[:6])",      TEAL),
        (0, 1, donor12,  "5' donor  12 bp (intron[:12])",     TEAL),
        (1, 0, accept6,  "3' acceptor  6 bp  (intron[-6:])",  RUST),
        (1, 1, accept12, "3' acceptor  12 bp (intron[-12:])", RUST),
    ]
    for r, c, series, title, color in cells:
        ax = axes[r, c]
        vc = series.value_counts().head(TOP_N)
        if vc.empty:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes, color="#999")
            ax.set_xticks([]); ax.set_yticks([])
            for sp_ in ax.spines.values(): sp_.set_visible(False)
        else:
            labels = vc.index.tolist()[::-1]
            vals = vc.values.tolist()[::-1]
            ax.barh(labels, vals, color=color)
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("intron count")
    out = os.path.join(outdir, f"{_slug(out_name or label)}.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--introns", required=True,
                    help="TSV with first40, last40, genome_id, organism")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label", default="introns",
                    help="label fragment in the suptitle")
    ap.add_argument("--min-records", type=int, default=5)
    ap.add_argument("--filter-col", default=None,
                    help="optional column to filter rows by before plotting")
    ap.add_argument("--filter-value", default=None,
                    help="rows must have filter-col == this value")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.introns, sep="\t", low_memory=False)
    df = df.dropna(subset=["first40", "last40"])
    df["first40"] = df["first40"].astype(str)
    df["last40"]  = df["last40"].astype(str)
    if args.filter_col and args.filter_value is not None:
        before = len(df)
        vals = [v for v in str(args.filter_value).split(",") if v]
        df = df[df[args.filter_col].astype(str).isin(vals)]
        print(f"filter {args.filter_col} in {vals}: {before} → {len(df)}",
              flush=True)

    plot_one("ALL_species", df, args.outdir,
             title_suffix=f"   {args.label}")
    print(f"wrote ALL_species  (n={len(df)})", flush=True)

    for (gid, organism), sub in df.groupby(["genome_id", "organism"], sort=False):
        if len(sub) < args.min_records:
            continue
        # Match the slug used by extract_telotron_fasta / splice-logo dirs:
        # "{gid}_{organism}".
        out_name = f"{_slug(gid)}_{_slug(organism)}"
        title = f"{organism} ({gid})"
        plot_one(out_name, sub, args.outdir, title_suffix=f"   {args.label}",
                 out_name=out_name, title=title)
        print(f"  {gid}: n={len(sub)}", flush=True)


if __name__ == "__main__":
    main()
