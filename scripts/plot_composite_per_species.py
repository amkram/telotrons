#!/usr/bin/env python3
"""Composite per-species figures: arbitrary N panels side by side.

Given N (dir, label) pairs, for each species filename (PNG) that appears in at
least one of the input dirs, compose a single horizontal figure containing one
panel per input. Missing PNGs in a given dir render as an empty "(no data)"
placeholder so categories without observed loci still show up explicitly.

The species universe is the union of filenames across all input dirs (a
species missing from one dir still appears, with that panel blank).
"""
import argparse
import os

import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def list_pngs(d):
    if not os.path.isdir(d):
        return set()
    return {f for f in os.listdir(d) if f.endswith(".png")}


def compose_one(paths, labels, out_path):
    n_panels = len(paths)
    # Read present panels; record dims so we can sensibly size missing ones.
    imgs = []
    for p in paths:
        imgs.append(mpimg.imread(p) if p and os.path.exists(p) else None)
    # Reference shape: max height among present panels (fallback 600).
    heights = [img.shape[0] for img in imgs if img is not None] or [600]
    target_h = max(heights)
    widths = []
    for img in imgs:
        if img is None:
            widths.append(target_h * 0.9)  # square-ish placeholder
        else:
            widths.append(img.shape[1] * (target_h / img.shape[0]))
    total_w = sum(widths)
    fig_w = total_w / 100
    fig_h = target_h / 100 + 0.4
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_w, fig_h),
                             layout="constrained",
                             gridspec_kw={"width_ratios": widths})
    if n_panels == 1:
        axes = [axes]
    for ax, img, label in zip(axes, imgs, labels):
        if img is None:
            ax.text(0.5, 0.5, "no loci", ha="center", va="center",
                    transform=ax.transAxes, color="#999", fontsize=12)
        else:
            ax.imshow(img)
        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.axis("off")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", action="append", required=True,
                    metavar="DIR:LABEL",
                    help="repeat once per panel; LABEL may contain colons after the first")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    panels = []
    for spec in args.panel:
        if ":" not in spec:
            raise SystemExit(f"--panel must be DIR:LABEL ({spec!r})")
        d, lbl = spec.split(":", 1)
        panels.append((d, lbl))

    os.makedirs(args.outdir, exist_ok=True)
    universe = set()
    per_dir = []
    for d, _ in panels:
        files = list_pngs(d)
        per_dir.append(files)
        universe |= files
    if not universe:
        print("no figures in any panel; nothing to compose")
        return
    print(f"composing {len(universe)} species across {len(panels)} panels",
          flush=True)

    for fname in sorted(universe):
        paths = [
            (os.path.join(d, fname) if fname in files else None)
            for (d, _), files in zip(panels, per_dir)
        ]
        labels = [lbl for _, lbl in panels]
        out = os.path.join(args.outdir, fname)
        compose_one(paths, labels, out)
        present = sum(1 for p in paths if p)
        print(f"  {fname}  ({present}/{len(panels)} panels present)", flush=True)


if __name__ == "__main__":
    main()
