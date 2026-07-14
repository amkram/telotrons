#!/usr/bin/env python3
"""Per-species figures: top boundary 6-mers (GTXXXX and XXXXAG) split by architecture.

Architectures plotted (rows): canonical set from `_common.ARCH_ORDER`
(includes `Other` so loci with no canonical architecture are not silently dropped).
Sides plotted (cols): donor (GT-XXXX) | linker_5' | linker_3' | acceptor (XXXX-AG).

NB: ranking is by **raw count** within (genome, architecture, side) — the input
`boundary_kmers_by_architecture.tsv` carries no enrichment statistic
(`analyze_telotrons.py` emits enrichment to a separate `boundary_kmer_enrichment.tsv`).
Plot is descriptive frequency only; titles reflect that.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

from _common import ARCH_COLORS, ARCH_ORDER, slug as _slug

TOP_N = 8


# Per-architecture active side set. `Other` defaults to donor/acceptor only;
# anything missing from this map also falls back to donor/acceptor.
SIDES_BY_ARCH = {
    "GT-F-AG":          ["donor", "acceptor"],
    "GT-R-AG":          ["donor", "acceptor"],
    "GT-F-R-AG":        ["donor", "acceptor"],
    "GT-F-linker-R-AG": ["donor", "linker_left", "linker_right", "acceptor"],
    "GT-R-linker-F-AG": ["donor", "linker_left", "linker_right", "acceptor"],
    "Other":            ["donor", "acceptor"],
}


def plot_one(genome_id, organism, sub, outdir):
    """sub: kmer table rows for this genome_id only."""
    n_rows = len(ARCH_ORDER)
    n_cols = 4  # donor | linker_left | linker_right | acceptor
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 1.6 * n_rows + 0.6),
                             layout="constrained")

    total_loci = sub[sub.side.isin(["donor"])]["count"].sum()
    # Descriptive frequency only — see module docstring for why this is not an enrichment plot.
    fig.suptitle(f"{organism}  ({genome_id})   {total_loci} loci   "
                 f"— most frequent boundary 6-mers (descriptive)",
                 fontsize=11, fontweight="bold")

    col_titles = ["donor  GT-XXXX", "linker_5'  XX-XXXX", "linker_3'  XXXX-XX", "acceptor  XXXX-AG"]
    side_for_col = ["donor", "linker_left", "linker_right", "acceptor"]

    for i, arch in enumerate(ARCH_ORDER):
        color = ARCH_COLORS.get(arch, "#BBBBBB")
        n_arch = sub[(sub.architecture == arch) & (sub.side == "donor")]["count"].sum()
        active_sides = SIDES_BY_ARCH.get(arch, ["donor", "acceptor"])
        for j in range(n_cols):
            ax = axes[i, j]
            side = side_for_col[j]
            if side not in active_sides:
                ax.text(0.5, 0.5, "n/a", ha="center", va="center",
                        color="#EEEEEE", transform=ax.transAxes, fontsize=10)
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                if j == 0:
                    ax.set_ylabel(f"{arch}\n(n={n_arch})", fontsize=9)
                if i == 0:
                    ax.set_title(col_titles[j], fontsize=10)
                continue
            d = sub[(sub.architecture == arch) & (sub.side == side)]
            if d.empty:
                ax.text(0.5, 0.5, "—", ha="center", va="center",
                        color="#BBBBBB", transform=ax.transAxes, fontsize=14)
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
            else:
                d = d.sort_values("count", ascending=False).head(TOP_N)
                labels = d.kmer.tolist()[::-1]
                vals = d["count"].tolist()[::-1]
                ax.barh(labels, vals, color=color)
                ax.tick_params(axis="y", labelsize=7)
                ax.tick_params(axis="x", labelsize=8)
                ax.spines[["top", "right"]].set_visible(False)
                ax.set_xlabel("count (descriptive)" if i == n_rows - 1 else "")
            if j == 0:
                ax.set_ylabel(f"{arch}\n(n={n_arch})", fontsize=9)
            if i == 0:
                ax.set_title(col_titles[j], fontsize=10)

    out = f"{outdir}/{_slug(genome_id)}_{_slug(organism)}.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmers", required=True, help="boundary_kmers_by_architecture.tsv")
    ap.add_argument("--species", required=True, help="final_species_summary.tsv (for organism names)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    km = pd.read_csv(args.kmers, sep="\t")
    sp = pd.read_csv(args.species, sep="\t")
    org_map = dict(zip(sp.genome_id, sp.organism))

    for gid, sub in km.groupby("genome_id"):
        organism = org_map.get(gid, gid)
        plot_one(gid, organism, sub, args.outdir)
        print(f"  {gid}")


if __name__ == "__main__":
    main()
