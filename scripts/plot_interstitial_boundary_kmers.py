#!/usr/bin/env python3
"""Per-species figures: top 5'- and 3'-end k-mers of interstitial telomeric arrays.

Layout per species (one PNG):
    2 rows (5' end, 3' end) × 2 cols (6-mer, 12-mer)
    Each panel: top-N k-mers ranked by count
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

from _common import slug as _slug, TEAL, RUST

TOP_N = 10


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
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.arrays, sep="\t", low_memory=False)
    if args.filter_col and args.filter_value is not None:
        before = len(df)
        vals = [v for v in str(args.filter_value).split(",") if v]
        df = df[df[args.filter_col].astype(str).isin(vals)]
        print(f"filter {args.filter_col} in {vals}: {before} → {len(df)}",
              flush=True)
    if df.empty:
        print("no interstitial arrays; nothing to plot")
        return
    for gid, sub in df.groupby("genome_id"):
        organism = sub.organism.iloc[0]
        plot_one(gid, organism, sub, args.outdir)
        print(f"  {gid}: {len(sub)}")
    print(f"wrote {df.genome_id.nunique()} figures into {args.outdir}/")


if __name__ == "__main__":
    main()
