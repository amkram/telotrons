#!/usr/bin/env python3
"""One figure per pipeline stage, summarising what that stage produced.

Stages plotted:
    1. manifests          -> stage1_manifests.png  (genomes per group)
    2. scan_all           -> stage3_scan.png       (introns scanned, raw candidates per species)
    3. filter_final       -> stage4_filter_funnel.png (raw -> filtered, per species)
    4. analyze            -> stage5_distance.png   (telotron vs control distance-to-end)
    5. architecture       -> stage6_architecture.png (per-species stacked architecture)
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import TEAL, RUST, GREY, ARCH_ORDER, ARCH_COLORS


def _label(row):
    return f"{row.organism}\n({row.genome_id})"


def plot_manifests(manifests, outpath):
    counts = manifests.groupby("group").size().sort_values(ascending=True)
    fig, ax = plt.subplots(figsize=(7, max(2.5, 0.5 * len(counts) + 1)), layout="constrained")
    ax.barh(counts.index, counts.values, color=TEAL)
    for y, v in enumerate(counts.values):
        ax.text(v + 0.05, y, str(v), va="center", fontsize=9)
    ax.set_xlabel("genomes")
    ax.set_title("Stage 1 — manifests: genomes per group")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_scan(raw, outpath):
    """Bar chart: introns scanned + raw candidates per species."""
    d = raw[raw.status == "OK"].copy()
    d = d.sort_values("telotron_candidates", ascending=True)
    labels = [_label(r) for _, r in d.iterrows()]
    y = np.arange(len(d))

    fig, axes = plt.subplots(1, 2, figsize=(12, max(3, 0.4 * len(d) + 1)),
                             layout="constrained", sharey=True)
    axes[0].barh(y, d.introns_scanned, color=TEAL)
    axes[0].set_yticks(y); axes[0].set_yticklabels(labels, fontsize=8)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("introns scanned (log)")
    axes[0].set_title("Stage 3 — scan: introns parsed by gt + scanned by seqkit")
    axes[0].spines[["top", "right"]].set_visible(False)

    axes[1].barh(y, d.telotron_candidates.clip(lower=0.5), color=RUST)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("raw candidates (log)")
    axes[1].set_title("Stage 3 — raw telotron candidates")
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_filter_funnel(raw, species_final, outpath):
    """For each species: raw candidates -> final filtered telotrons."""
    raw_ok = raw[raw.status == "OK"][["genome_id", "organism", "telotron_candidates"]]
    fin = species_final[["genome_id", "telotrons"]]
    m = raw_ok.merge(fin, on="genome_id", how="left").fillna({"telotrons": 0})
    m = m[m.telotron_candidates > 0].copy()
    m["retention"] = m.telotrons / m.telotron_candidates.replace(0, np.nan)
    m = m.sort_values("telotrons", ascending=True)
    labels = [_label(r) for _, r in m.iterrows()]
    y = np.arange(len(m))

    fig, ax = plt.subplots(figsize=(11, max(3, 0.45 * len(m) + 1)), layout="constrained")
    ax.barh(y, m.telotron_candidates, color=GREY, label="raw candidates")
    ax.barh(y, m.telotrons, color=TEAL, label="final telotrons")
    for yi, (rc, ft, ret) in enumerate(zip(m.telotron_candidates, m.telotrons, m.retention)):
        pct = f"{100 * ret:.1f}%" if pd.notna(ret) else "—"
        ax.text(max(rc, 1) * 1.1, yi, f"{int(ft)} / {int(rc)}  ({pct})",
                va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("count (log)")
    ax.set_title("Stage 4 — filter: raw → filtered telotrons per species")
    ax.legend(loc="lower right", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_distance(dist, raw, outpath):
    """Telotron median distance-to-end vs length-matched control, per species."""
    d = dist.merge(raw[["genome_id", "organism"]], on="genome_id", how="left")
    d = d[d.telotron_n > 0].copy()
    d = d.sort_values("median_telotron_distance_to_end")
    y = np.arange(len(d))
    labels = [_label(r) for _, r in d.iterrows()]

    fig, ax = plt.subplots(figsize=(10, max(3, 0.45 * len(d) + 1)), layout="constrained")
    ax.scatter(d.median_telotron_distance_to_end / 1e3, y, color=TEAL, label="telotron median", s=40)
    ax.scatter(d.median_control_distance_to_end / 1e3, y, color=GREY, label="control median", s=40)
    for yi, (t, c) in enumerate(zip(d.median_telotron_distance_to_end, d.median_control_distance_to_end)):
        ax.plot([t / 1e3, c / 1e3], [yi, yi], color="#CCCCCC", linewidth=1, zorder=0)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("distance to contig end (kb, log)")
    ax.set_title("Stage 5 — analyze: telotron vs length-matched control distance-to-end")
    ax.legend(loc="lower right", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def plot_architecture(arch_df, outpath):
    """Stacked bar per species: counts in each of the 5 architecture categories."""
    by = arch_df.groupby(["genome_id", "organism", "architecture"]).size().unstack(fill_value=0)
    for a in ARCH_ORDER:
        if a not in by.columns:
            by[a] = 0
    by = by[ARCH_ORDER]
    by["_total"] = by.sum(axis=1)
    by = by.sort_values("_total", ascending=True)
    by = by.drop(columns=["_total"])
    labels = [f"{org}\n({gid})" for gid, org in by.index]
    y = np.arange(len(by))

    fig, ax = plt.subplots(figsize=(11, max(3, 0.5 * len(by) + 1)), layout="constrained")
    left = np.zeros(len(by))
    for a in ARCH_ORDER:
        vals = by[a].values
        ax.barh(y, vals, left=left, color=ARCH_COLORS[a], label=a, edgecolor="white", linewidth=0.3)
        left += vals
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("telotron loci (log)")
    ax.set_title("Stage 6 — architecture: F/R/linker composition per species")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifests", required=True)
    ap.add_argument("--raw-summary", required=True)
    ap.add_argument("--species-final", required=True)
    ap.add_argument("--distance", required=True)
    ap.add_argument("--architecture-loci", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    manifests = pd.read_csv(args.manifests, sep="\t")
    raw = pd.read_csv(args.raw_summary, sep="\t")
    species_final = pd.read_csv(args.species_final, sep="\t")
    dist = pd.read_csv(args.distance, sep="\t")
    arch = pd.read_csv(args.architecture_loci, sep="\t")

    plot_manifests(manifests, f"{args.outdir}/stage1_manifests.png")
    plot_scan(raw, f"{args.outdir}/stage3_scan.png")
    plot_filter_funnel(raw, species_final, f"{args.outdir}/stage4_filter_funnel.png")
    plot_distance(dist, raw, f"{args.outdir}/stage5_distance.png")
    plot_architecture(arch[arch.architecture != "Unknown"], f"{args.outdir}/stage6_architecture.png")
    print(f"wrote 5 stage figures into {args.outdir}/")


if __name__ == "__main__":
    main()
