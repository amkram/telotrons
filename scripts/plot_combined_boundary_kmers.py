#!/usr/bin/env python3
"""Compose a single figure with intron + interstitial boundary k-mer panels.

For every species (taken from interstitial directory listing or arch directory
listing), pastes its intron-architecture boundary panel on the left and its
interstitial boundary panel on the right. Rows ordered by `group` (mammals,
plants, etc.) then by organism name so related taxa cluster.

Output: a multi-page PDF (rows split into pages so each row is legible) and a
single tall PNG ("--png-too").
"""
import argparse
import os

import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from _common import slug as _slug


def collect(species_tsv, intron_dir, interst_dir):
    sp_df = pd.read_csv(species_tsv, sep="\t")
    rows = []
    for r in sp_df.itertuples(index=False):
        slug = f"{_slug(r.genome_id)}_{_slug(r.organism)}"
        intr = os.path.join(intron_dir, f"{slug}.png")
        inte = os.path.join(interst_dir, f"{slug}.png")
        if not (os.path.exists(intr) or os.path.exists(inte)):
            continue
        rows.append({
            "genome_id": r.genome_id, "organism": r.organism,
            "group": getattr(r, "group", "") if hasattr(r, "group") else "",
            "intron_png": intr if os.path.exists(intr) else None,
            "interst_png": inte if os.path.exists(inte) else None,
        })
    rows.sort(key=lambda x: (str(x.get("group", "")), x["organism"]))
    return rows


def make_pdf(rows, out_pdf, rows_per_page=2):
    """Each page: `rows_per_page` species. Two columns: intron | interstitial."""
    with PdfPages(out_pdf) as pdf:
        for i in range(0, len(rows), rows_per_page):
            chunk = rows[i:i + rows_per_page]
            fig, axes = plt.subplots(
                len(chunk), 2,
                figsize=(22, 7.5 * len(chunk)),
                layout="constrained",
                squeeze=False,
            )
            for r_idx, rec in enumerate(chunk):
                title = f"{rec['organism']}  ({rec['genome_id']})   group={rec.get('group') or '—'}"
                for c_idx, (key, label) in enumerate(
                    [("intron_png", "intron telotron boundary k-mers (by architecture)"),
                     ("interst_png", "interstitial array boundary k-mers (non-genic, non-intronic)")]
                ):
                    ax = axes[r_idx, c_idx]
                    p = rec[key]
                    if p and os.path.exists(p):
                        img = Image.open(p)
                        ax.imshow(img)
                    else:
                        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                                color="#aaa", transform=ax.transAxes, fontsize=14)
                    ax.set_xticks([]); ax.set_yticks([])
                    for sp_ in ax.spines.values(): sp_.set_visible(False)
                    if r_idx == 0:
                        ax.set_title(label, fontsize=12, fontweight="bold")
                # vertical species label on the leftmost
                axes[r_idx, 0].set_ylabel(title, fontsize=11, fontweight="bold",
                                           rotation=90, labelpad=10)
            pdf.savefig(fig, dpi=140)
            plt.close(fig)


def make_tall_png(rows, out_png):
    """Concatenate every species row into one tall PNG (intron | interstitial)."""
    pairs = []
    for rec in rows:
        a = Image.open(rec["intron_png"]) if rec["intron_png"] else None
        b = Image.open(rec["interst_png"]) if rec["interst_png"] else None
        if a is None and b is None:
            continue
        # normalise heights
        h = max((img.height for img in (a, b) if img is not None))
        if a is not None and a.height != h:
            a = a.resize((int(a.width * h / a.height), h))
        if b is not None and b.height != h:
            b = b.resize((int(b.width * h / b.height), h))
        # fill missing with blank
        if a is None: a = Image.new("RGB", (b.width, h), "white")
        if b is None: b = Image.new("RGB", (a.width, h), "white")
        row = Image.new("RGB", (a.width + b.width + 20, h), "white")
        row.paste(a, (0, 0))
        row.paste(b, (a.width + 20, 0))
        pairs.append((rec, row))

    # Stack vertically with a small label band per row
    if not pairs:
        return
    max_w = max(r.width for _, r in pairs)
    band_h = 36
    total_h = sum(r.height + band_h for _, r in pairs)
    out = Image.new("RGB", (max_w, total_h), "white")
    y = 0
    try:
        from PIL import ImageDraw, ImageFont
        font = None
        for path in ("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                     "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf"):
            if os.path.exists(path):
                font = ImageFont.truetype(path, 22); break
        if font is None: font = ImageFont.load_default()
    except Exception:
        font, ImageDraw = None, None

    for rec, row in pairs:
        if ImageDraw is not None:
            band = Image.new("RGB", (max_w, band_h), (235, 235, 235))
            d = ImageDraw.Draw(band)
            label = f"{rec['organism']}  ({rec['genome_id']})   group={rec.get('group') or '—'}"
            d.text((10, 6), label, fill="black", font=font)
            out.paste(band, (0, y))
        y += band_h
        out.paste(row, (0, y))
        y += row.height
    out.save(out_png, optimize=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True, help="final_species_summary.tsv")
    ap.add_argument("--intron-dir", required=True,
                    help="figures/boundary_kmers_by_arch directory")
    ap.add_argument("--interstitial-dir", required=True,
                    help="figures/interstitial_boundary_kmers directory")
    ap.add_argument("--out-pdf", required=True)
    ap.add_argument("--out-png", required=False, default=None)
    ap.add_argument("--rows-per-page", type=int, default=2)
    args = ap.parse_args()

    rows = collect(args.species, args.intron_dir, args.interstitial_dir)
    if not rows:
        print("no species rows; nothing to compose"); return
    make_pdf(rows, args.out_pdf, args.rows_per_page)
    print(f"wrote {args.out_pdf}  ({len(rows)} species, {args.rows_per_page} per page)")
    if args.out_png:
        make_tall_png(rows, args.out_png)
        print(f"wrote {args.out_png}")


if __name__ == "__main__":
    main()
