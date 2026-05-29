#!/usr/bin/env python3
"""Sequence logos of interstitial telomeric array boundaries.

For each species (and an overall combined panel), build two logos:
  5' end: 10 bp upstream flank ++ first repeat unit
  3' end: last repeat unit     ++ 10 bp downstream flank

Sequences are in the canonical (G-rich) telomeric orientation as recorded by
find_interstitial_arrays.py. Because units come in three lengths (TTAGG=5,
TTAGGG=6, TTTAGGG=7), per-species logos pad the unit side with 'N' so columns
align across arrays; the unit block is right-padded on the 5' panel and
left-padded on the 3' panel. Logos use logomaker information-bit heights.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import logomaker

from _common import slug as _slug

ACGTN = list("ACGT")


def stack(seqs, length, side):
    """Return a position-frequency DataFrame of width `length`.

    side='right' pads short strings with 'N' on the right (used for the 3'
    unit side on the 5' panel); side='left' pads on the left.
    Ns are ignored when building the count matrix.
    """
    counts = {b: [0] * length for b in ACGTN}
    for s in seqs:
        s = s.upper()
        if side == "right":
            s = s + "N" * (length - len(s))
        else:
            s = "N" * (length - len(s)) + s
        s = s[:length] if side == "right" else s[-length:]
        for i, ch in enumerate(s):
            if ch in counts:
                counts[ch][i] += 1
    df = pd.DataFrame(counts)
    # normalise rows; if a column is all-N, leave it as zeros
    row_sums = df.sum(axis=1).replace(0, 1)
    return df.div(row_sums, axis=0)


def to_information(freq_df):
    """Convert per-position frequencies to information-bit heights."""
    return logomaker.transform_matrix(
        freq_df, from_type="probability", to_type="information"
    )


def draw_panel(ax, seqs, side, flank_len, unit_len, title, color_scheme="classic"):
    if side == "5p":
        # upstream flank (left, fixed=flank_len) ++ first unit (right, variable)
        total = flank_len + unit_len
        flank_part = [s[:flank_len] for s in seqs]
        unit_part  = [s[flank_len:flank_len + unit_len] for s in seqs]
        flank_df = stack(flank_part, flank_len, "left")
        unit_df  = stack(unit_part,  unit_len, "right")
        full_df  = pd.concat([flank_df, unit_df], ignore_index=True)
    else:
        # last unit (left, variable) ++ downstream flank (right, fixed)
        total = unit_len + flank_len
        unit_part  = [s[:unit_len] for s in seqs]
        flank_part = [s[unit_len:unit_len + flank_len] for s in seqs]
        unit_df  = stack(unit_part,  unit_len, "left")
        flank_df = stack(flank_part, flank_len, "right")
        full_df  = pd.concat([unit_df, flank_df], ignore_index=True)
    info = to_information(full_df)
    logomaker.Logo(info, ax=ax, color_scheme=color_scheme,
                   show_spines=False, fade_below=0.0)
    ax.set_ylim(0, 2)
    ax.set_yticks([0, 1, 2])
    ax.set_ylabel("bits")
    ax.set_title(title, fontsize=10)
    # mark unit/flank boundary
    if side == "5p":
        ax.axvline(flank_len - 0.5, color="#888", lw=0.8, linestyle="--")
    else:
        ax.axvline(unit_len - 0.5, color="#888", lw=0.8, linestyle="--")
    ax.set_xticks(range(total))
    ax.set_xticklabels([str(i) for i in range(total)], fontsize=7)


def build_logo_inputs(sub):
    """Return (five_prime_strings, three_prime_strings, unit_len)."""
    sub = sub[(sub.first_unit != "") & (sub.last_unit != "")]
    if sub.empty:
        return [], [], 0
    unit_len = int(sub.first_unit.str.len().mode().iat[0])
    sub = sub[sub.first_unit.str.len() == unit_len]
    sub = sub[sub.last_unit.str.len() == unit_len]
    five = (sub.upstream_flank10 + sub.first_unit).tolist()
    three = (sub.last_unit + sub.downstream_flank10).tolist()
    return five, three, unit_len


def plot_one(label, sub, outdir, flank_len=10, out_name=None):
    five, three, unit_len = build_logo_inputs(sub)
    n = len(five)
    if n == 0 or unit_len == 0:
        return False
    fig, axes = plt.subplots(2, 1, figsize=(0.45 * (flank_len + unit_len) + 2.5, 4.6),
                             layout="constrained")
    fig.suptitle(f"{label}   n={n}   unit={unit_len} bp",
                 fontsize=11, fontweight="bold")
    draw_panel(axes[0], five, "5p", flank_len, unit_len,
               "5' end: upstream flank (left of |) + first repeat unit")
    draw_panel(axes[1], three, "3p", flank_len, unit_len,
               "3' end: last repeat unit + downstream flank (right of |)")
    out = os.path.join(outdir, f"{_slug(out_name or label)}.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arrays", required=True, help="results/interstitial_arrays.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank-len", type=int, default=10)
    ap.add_argument("--min-arrays", type=int, default=5,
                    help="skip species with fewer arrays than this")
    ap.add_argument("--filter-col", default=None)
    ap.add_argument("--filter-value", default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.arrays, sep="\t")
    needed = {"upstream_flank10", "first_unit", "last_unit", "downstream_flank10"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"interstitial_arrays.tsv missing columns: {missing}")
    df = df.fillna({"upstream_flank10": "", "first_unit": "",
                    "last_unit": "", "downstream_flank10": ""})
    if args.filter_col and args.filter_value is not None:
        before = len(df)
        vals = [v for v in str(args.filter_value).split(",") if v]
        df = df[df[args.filter_col].astype(str).isin(vals)]
        print(f"filter {args.filter_col} in {vals}: {before} → {len(df)}",
              flush=True)

    # overall (all species)
    plot_one("ALL_species", df, args.outdir, flank_len=args.flank_len)
    print(f"wrote ALL_species logo  (n={len(df)})", flush=True)

    # per species
    for (gid, organism), sub in df.groupby(["genome_id", "organism"], sort=False):
        if len(sub) < args.min_arrays:
            continue
        label = f"{organism} ({gid})"
        out_name = f"{_slug(gid)}_{_slug(organism)}"
        ok = plot_one(label, sub, args.outdir, flank_len=args.flank_len,
                      out_name=out_name)
        if ok:
            print(f"  {gid}: n={len(sub)}", flush=True)


if __name__ == "__main__":
    main()
