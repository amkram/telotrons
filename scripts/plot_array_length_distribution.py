#!/usr/bin/env python3
"""Per-species length distributions: interstitial telomeric arrays vs telotrons.

One row per species (those with data in at least one input); two columns
(ITS array_len, telotron telomeric_bases). Single PNG."""
import argparse
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_by_species(path, col, min_len=0):
    out = defaultdict(list)
    names = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            try:
                v = int(r[col])
            except (KeyError, ValueError):
                continue
            if v < max(1, min_len):
                continue
            gid = r.get("genome_id", "")
            out[gid].append(v)
            names.setdefault(gid, r.get("organism", gid))
    return {k: np.array(v) for k, v in out.items()}, names


# Shared bin edges: 25-bp bins 0-1000, plus one overflow bin for >1000 bp.
BIN_WIDTH = 25
BIN_EDGES = list(range(0, 1001, BIN_WIDTH)) + [float("inf")]
N_BINS = len(BIN_EDGES) - 1
BIN_CENTERS = [(BIN_EDGES[i] + BIN_WIDTH / 2) for i in range(N_BINS - 1)] + [1000 + BIN_WIDTH / 2]


def panel(ax, vals, color, ylabel=None):
    if len(vals) == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        return
    counts, _ = np.histogram(vals, bins=BIN_EDGES)
    ax.bar(BIN_CENTERS, counts, width=BIN_WIDTH, color=color,
           edgecolor="black", linewidth=0.2, align="center")
    ax.set_xlim(0, 1000 + BIN_WIDTH)
    ax.set_xticks([0, 200, 400, 600, 800, 1000, 1000 + BIN_WIDTH / 2])
    ax.set_xticklabels(["0", "200", "400", "600", "800", "1000", ">1k"])
    med = int(np.median(vals))
    p95 = int(np.percentile(vals, 95))
    ax.set_title(f"n={len(vals)}  med={med}  p95={p95}",
                 fontsize=8, loc="left")
    ax.tick_params(labelsize=7)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interstitial", required=True)
    ap.add_argument("--telotrons", required=True)
    ap.add_argument("--non-telotrons", required=True,
                    help="results/non_telotron_controls.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=0,
                    help="drop arrays/telotrons shorter than this (bp)")
    args = ap.parse_args()

    its, its_names = read_by_species(args.interstitial, "array_len",
                                      min_len=args.min_len)
    telo, telo_names = read_by_species(args.telotrons, "telomeric_bases",
                                        min_len=args.min_len)
    # Non-telotron introns: length = intron_len (no telomeric array exists).
    # min_len filter still applied for consistency.
    nt, nt_names = read_by_species(args.non_telotrons, "intron_len",
                                    min_len=args.min_len)

    def keep(gid):
        return gid.startswith("GCF_000499") or gid.startswith("TARA_")

    species = sorted([g for g in (set(its) | set(telo) | set(nt)) if keep(g)],
                     key=lambda g: -(len(its.get(g, [])) + len(telo.get(g, []))
                                       + len(nt.get(g, []))))
    if not species:
        raise SystemExit("no species with data")
    names = {**nt_names, **telo_names, **its_names}

    n = len(species)
    fig, axes = plt.subplots(n, 3, figsize=(13, 1.6 * n + 0.5), squeeze=False)
    for i, gid in enumerate(species):
        label = f"{names.get(gid, gid)}\n{gid}"
        panel(axes[i, 0], its.get(gid, np.array([])), "#4878a8", ylabel=label)
        panel(axes[i, 1], telo.get(gid, np.array([])), "#c45a3a")
        panel(axes[i, 2], nt.get(gid, np.array([])), "#6aa84f")
    axes[0, 0].set_title("Interstitial telomeric arrays (array_len)",
                          fontsize=10, loc="center", pad=14)
    axes[0, 1].set_title("Telotrons (telomeric_bases in intron)",
                          fontsize=10, loc="center", pad=14)
    axes[0, 2].set_title("Non-telotron introns (intron_len)",
                          fontsize=10, loc="center", pad=14)
    for c in range(3):
        axes[-1, c].set_xlabel("length (bp)", fontsize=9)

    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"wrote {args.out}  ({n} species)", flush=True)


if __name__ == "__main__":
    main()
