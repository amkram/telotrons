#!/usr/bin/env python3
"""Per-species telotron figures + one aggregate counts chart."""
import argparse
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

from _common import (
    slug as _slug,
    TEAL, RUST, GOLD, GREY,
    ORIENT_COLORS,
    ARCH_ORDER, ARCH_COLORS,
)


def plot_counts(species, outdir):
    """Aggregate: horizontal bar of telotron counts, all positive species."""
    d = species[species.telotrons > 0].sort_values("telotrons")
    if d.empty:
        return
    labels = [f"{r.organism}\n({r.genome_id})" for _, r in d.iterrows()]
    fig, ax = plt.subplots(figsize=(9, max(3, len(d) * 0.45 + 1)), layout="constrained")
    ax.barh(labels, d.telotrons.values, color=TEAL)
    ax.set_xscale("log")
    ax.set_xlabel("final telotron loci")
    ax.set_title("Telotron-positive genomes")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(f"{outdir}/telotron_counts.png", dpi=200)
    plt.close(fig)


def plot_species(genome_id, organism, final_g, kmers_g, dist_row, arch_g, outdir):
    """One 3-panel figure for a single positive species."""
    fig = plt.figure(figsize=(11, 7), layout="constrained")
    gs = gridspec.GridSpec(2, 2, figure=fig)

    # ── Panel A: splice architecture breakdown ──────────────────────────────
    # Drives the panel from `architecture` (classify_telotron_architecture)
    # rather than `orientation` so the bidirectional GT-F-R-AG /
    # GT-{F,R}-linker-{R,F}-AG classes are surfaced individually instead of
    # being collapsed into the orientation `Mixed` bin. See
    # verify_figures_visualization new_concern 4 (2026-06-04) +
    # telotron_two_linker_classes_2026-06-03. Falls back to `orientation` if
    # the `architecture` column is missing from older final-set TSVs.
    ax_or = fig.add_subplot(gs[0, 0])
    if "architecture" in final_g.columns and final_g["architecture"].notna().any():
        arch_counts = final_g["architecture"].value_counts()
        ordered = [a for a in ARCH_ORDER if a in arch_counts.index]
        for a in arch_counts.index:
            if a not in ordered:
                ordered.append(a)
        arch_counts = arch_counts.reindex(ordered)
        colors = [ARCH_COLORS.get(a, GREY) for a in arch_counts.index]
        ax_or.bar(arch_counts.index, arch_counts.values, color=colors)
        ax_or.set_title("Splice architecture")
    else:
        orient_counts = final_g["orientation"].value_counts()
        colors = [ORIENT_COLORS.get(o, GREY) for o in orient_counts.index]
        ax_or.bar(orient_counts.index, orient_counts.values, color=colors)
        ax_or.set_title("Splice orientation")
    ax_or.set_ylabel("loci")
    ax_or.spines[["top", "right"]].set_visible(False)
    for tick in ax_or.get_xticklabels():
        tick.set_rotation(30)
        tick.set_ha("right")

    # ── Panel B: intron length distribution ─────────────────────────────────
    ax_len = fig.add_subplot(gs[0, 1])
    lengths = final_g["intron_len"].dropna()
    if len(lengths):
        ax_len.hist(lengths, bins=min(40, max(10, len(lengths) // 5)),
                    color=TEAL, edgecolor="white", linewidth=0.3)
    ax_len.set_xlabel("intron length (bp)")
    ax_len.set_ylabel("loci")
    ax_len.set_title("Intron length")
    ax_len.spines[["top", "right"]].set_visible(False)

    # ── Panel C: top boundary k-mers ────────────────────────────────────────
    # analyze_telotrons.py writes Fisher p + BH q (`bh_q`) and flags the
    # boundary-k-mer comparison as "partly DEFINITIONAL" (telomere-motif
    # rotations are guaranteed enriched). Filter by q<0.05 when the column
    # is present; otherwise retitle as descriptive so the reader cannot
    # misread the panel as an enrichment test
    # (verify_figures_visualization Finding 7, 2026-06-04).
    ax_km = fig.add_subplot(gs[1, 0])
    if kmers_g is not None and not kmers_g.empty:
        if "bh_q" in kmers_g.columns:
            sig = kmers_g[kmers_g["bh_q"] < 0.05]
            title = "Boundary k-mers (q<0.05, top 15)"
        else:
            sig = kmers_g
            title = "Boundary k-mers (top 15, descriptive)"
        if not sig.empty:
            top = sig.sort_values("fold_enrichment", ascending=False).head(15)
            ax_km.barh(top["kmer"][::-1], top["fold_enrichment"][::-1], color=RUST)
            ax_km.set_xlabel("fold enrichment")
            n_shown = len(top)
            if n_shown < 15:
                title = f"{title.split(',')[0]} (n={n_shown})"
            ax_km.set_title(title)
            ax_km.spines[["top", "right"]].set_visible(False)
        else:
            ax_km.text(0.5, 0.5, "no kmers pass q<0.05", ha="center", va="center",
                       transform=ax_km.transAxes, color=GREY)
            ax_km.set_title("Boundary k-mers")
            ax_km.axis("off")
    else:
        ax_km.text(0.5, 0.5, "no boundary k-mer data", ha="center", va="center",
                   transform=ax_km.transAxes, color=GREY)
        ax_km.set_title("Boundary k-mers")
        ax_km.axis("off")

    # ── Panel D: distance-to-contig-end histogram ───────────────────────────
    ax_dist = fig.add_subplot(gs[1, 1])
    dists = final_g["distance_to_end"].dropna()
    if len(dists):
        ax_dist.hist(dists / 1e3, bins=min(40, max(10, len(dists) // 5)),
                     color=GOLD, edgecolor="white", linewidth=0.3)
    ax_dist.set_xlabel("distance to contig end (kb)")
    ax_dist.set_ylabel("loci")
    ax_dist.set_title("Distance to contig end")
    ax_dist.spines[["top", "right"]].set_visible(False)
    if dist_row is not None:
        # The histogram is over every positive locus (all contigs); the
        # median line comes from analyze_telotrons.distance_to_end() which
        # restricts to capped contigs (capped_only=True). Label the line
        # so the two denominators are not visually conflated
        # (verify_figures_visualization Finding 6, 2026-06-04).
        med = dist_row["median_telotron_distance_to_end"] / 1e3
        ax_dist.axvline(med, color=RUST, linestyle="--", linewidth=1,
                        label=f"median {med:.1f} kb (capped contigs)")
        ax_dist.legend(fontsize=8)

    n = len(final_g)
    # Summarise all distinct terminal motifs at the species level rather
    # than reading `.iloc[0]` as if it were scalar — 18% of focal Eimeria
    # telotrons mix TTAGGG + TTTAGGG at the same locus
    # (telotron_within_locus_motif_mixing_2026-06-03;
    # verify_figures_visualization new_concern 5, 2026-06-04).
    if "terminal_motif" in final_g.columns:
        motifs = (
            final_g["terminal_motif"].dropna().astype(str)
            .replace({"": pd.NA}).dropna().unique().tolist()
        )
        motif = ",".join(sorted(motifs)) if motifs else ""
    else:
        motif = ""
    fig.suptitle(f"{organism}  ({genome_id})   n={n} loci   motif: {motif}",
                 fontsize=10, fontweight="bold")

    slug = _slug(genome_id)
    fig.savefig(f"{outdir}/per_species/{slug}_{_slug(organism)}.png", dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--final", required=True)
    ap.add_argument("--kmers", required=True)
    ap.add_argument("--distance", required=True)
    ap.add_argument("--architecture", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(f"{args.outdir}/per_species", exist_ok=True)

    species = pd.read_csv(args.species, sep="\t")
    final = pd.read_csv(args.final, sep="\t")
    kmers = pd.read_csv(args.kmers, sep="\t")
    dist = pd.read_csv(args.distance, sep="\t")

    # aggregate counts figure
    plot_counts(species, args.outdir)

    # per-species figures for every positive genome
    positive = species[species.telotrons > 0]
    for _, row in positive.iterrows():
        gid = row.genome_id
        final_g = final[final.genome_id == gid]
        kmers_g = kmers[kmers.genome_id == gid] if "genome_id" in kmers.columns else pd.DataFrame()
        dist_rows = dist[dist.genome_id == gid]
        dist_row = dist_rows.iloc[0] if len(dist_rows) else None
        arch_g = pd.DataFrame()
        plot_species(gid, row.organism, final_g, kmers_g, dist_row, arch_g, args.outdir)


if __name__ == "__main__":
    main()
