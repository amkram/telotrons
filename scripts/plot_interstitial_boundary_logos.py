#!/usr/bin/env python3
"""Sequence logos of interstitial telomeric array boundaries.

For each species (and an overall combined panel), build two logos:
  5' end: 10 bp upstream flank  (telomere-rotation MASKED) ++ first repeat unit
  3' end: last repeat unit                                  ++ 10 bp downstream flank (masked)

Sequences are in the canonical (G-rich) telomeric orientation as recorded by
find_interstitial_arrays.py. Multi-motif arrays (mixed unit lengths) are NOT
silently dropped — they are routed to per-unit-length panels (one logo per
length) so the TTAGGG+TTTAGGG majority on Eimeria stays visible
(telotron_within_locus_motif_mixing_2026-06-03).

Statistical correctness:
- Rows are first filtered to ITS gold standard (Giulotto: >=4 units, identity
  >= 6/7 to perfect tandem) via `tandem_identity()` to drop the 73% of raw
  interstitial_arrays.tsv rows that are degenerate microsatellite noise
  (its_detection_gold_standard).
- Flanks are masked with `mask_telomere_fragments()` (rotations + revcomps,
  >= 2 tandem units) BEFORE concatenation, so a telomere rotation bleeding
  into the flank does not appear as a flank consensus
  (eimeria_subtelomere_rlfs_2026-05-15 retraction).
- Per-genome ACGT background is used as the logomaker `background=` so that
  AT-rich genomes (Eimeria ~80% AT) do not get artificially inflated bits on
  every A/T column (Schneider & Stephens 1990 NAR 18:6097).
- Counts -> information directly via `from_type='counts'` so the Schneider
  1986 small-sample correction `e(n)` is applied (WebLogo default).
- Logos with n < MIN_N (default 20) are skipped (WebLogo's threshold).
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import logomaker

from _common import slug as _slug, expand_motifs, read_fasta, find_genome_fasta
from telomere_mask import mask_telomere_fragments

ACGT = list("ACGT")
ACGTN = list("ACGTN")

# WebLogo-style minimum sample threshold for refusing to draw a logo.
MIN_N_FOR_LOGO = 20

# Telomere motifs used for both flank-masking and gold-standard tandem identity.
TELO_MOTIFS = ("TTAGGG", "TTTAGGG", "TTAGG", "TTGGGG")

# Backstop uniform background when per-genome ACGT cannot be computed.
UNIFORM_BG = {b: 0.25 for b in ACGT}


def tandem_identity(seq, motif):
    """Best identity of `seq` to a PERFECT tandem `(motif)n` over all phases
    x {fwd, revcomp}. Giulotto gold standard: ITS requires len >= 4 units AND
    identity >= 6/7 (= <1 mismatch per unit). Rejects degenerate microsatellite
    noise that find_arrays' fuzzy <=1mm-per-unit matching admits.
    """
    s = seq.upper()
    L = len(s)
    if L < len(motif):
        return 0.0
    rc = motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    best = 0
    for base in (motif, rc):
        for ph in range(len(base)):
            rot = base[ph:] + base[:ph]
            ref = (rot * (L // len(rot) + 2))[:L]
            m = sum(1 for a, b in zip(s, ref) if a == b)
            if m > best:
                best = m
    return best / L


def apply_gold_standard(df, min_units=4, min_ident=6.0 / 7.0):
    """Drop rows whose called array fails the Giulotto >=4-unit / >=6/7-identity
    test. Uses `array_seq` if present, otherwise reconstructs from
    `first_unit` * `n_units` -ish. Returns the filtered frame and a count tuple.
    """
    if df.empty:
        return df, (0, 0)
    n_before = len(df)
    # Best-available sequence to test: prefer a stored array_seq, else fall back
    # to the first_unit field repeated (rough but always present) for length
    # check + identity check.
    if "array_seq" in df.columns:
        seqs = df["array_seq"].fillna("").astype(str)
    else:
        # Reconstruct a minimum-length proxy: assume n_units * len(first_unit).
        if "n_units" in df.columns:
            seqs = (df["first_unit"].fillna("").astype(str)
                    * df["n_units"].fillna(0).astype(int))
        else:
            seqs = df["first_unit"].fillna("").astype(str)
    keep = []
    for s, motif in zip(seqs, df["first_unit"].fillna("").astype(str)):
        if not motif:
            keep.append(False)
            continue
        if len(s) < min_units * len(motif):
            keep.append(False)
            continue
        keep.append(tandem_identity(s, motif) >= min_ident)
    out = df[pd.Series(keep, index=df.index)].copy()
    return out, (n_before, len(out))


def compute_acgt_background(fasta_path):
    """Return per-base frequency dict over ACGT for the assembly FASTA."""
    counts = {b: 0 for b in ACGT}
    seqs, _ = (read_fasta(fasta_path), None) if False else (None, None)
    # Stream FASTA to avoid loading whole assembly: read in chunks.
    from _common import open_maybe_gz
    with open_maybe_gz(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            for ch in line.strip().upper():
                if ch in counts:
                    counts[ch] += 1
    total = sum(counts.values())
    if total == 0:
        return dict(UNIFORM_BG)
    return {b: counts[b] / total for b in ACGT}


def stack_counts(seqs, length, side):
    """Return a position-COUNT DataFrame of width `length` over ACGT only.

    side='right' pads short strings with 'N' on the right; side='left' pads on
    the left. Ns are not counted toward any base, so they reduce per-column
    total counts (which feeds the Schneider 1986 e(n) small-sample correction
    when logomaker transforms counts->information).
    """
    counts = {b: [0] * length for b in ACGT}
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
    return pd.DataFrame(counts)


def counts_to_information(counts_df, background):
    """Counts -> information with the Schneider 1986 e(n) small-sample
    correction (applied by logomaker when `from_type='counts'`) and an explicit
    non-uniform background.
    """
    bg = [background[b] for b in ACGT]
    return logomaker.transform_matrix(
        counts_df,
        from_type="counts",
        to_type="information",
        background=bg,
        pseudocount=1,
    )


def draw_panel(ax, seqs, side, flank_len, unit_len, title, background,
               color_scheme="classic"):
    if side == "5p":
        # upstream flank (left, fixed=flank_len) ++ first unit (right, variable)
        total = flank_len + unit_len
        flank_part = [s[:flank_len] for s in seqs]
        unit_part = [s[flank_len:flank_len + unit_len] for s in seqs]
        flank_counts = stack_counts(flank_part, flank_len, "left")
        unit_counts = stack_counts(unit_part, unit_len, "right")
        full_counts = pd.concat([flank_counts, unit_counts], ignore_index=True)
    else:
        # last unit (left, variable) ++ downstream flank (right, fixed)
        total = unit_len + flank_len
        unit_part = [s[:unit_len] for s in seqs]
        flank_part = [s[unit_len:unit_len + flank_len] for s in seqs]
        unit_counts = stack_counts(unit_part, unit_len, "left")
        flank_counts = stack_counts(flank_part, flank_len, "right")
        full_counts = pd.concat([unit_counts, flank_counts], ignore_index=True)

    info = counts_to_information(full_counts, background)
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


def mask_flanks(df, flank_len):
    """Return df copy where upstream_flank10 / downstream_flank10 have telomere
    rotations + RCs + tandem >= 2 unit runs masked to N. Prevents telomere
    bleed-into-flank from appearing as a flank consensus
    (eimeria_subtelomere_rlfs_2026-05-15 retraction).
    """
    out = df.copy()
    out["upstream_flank10"] = out["upstream_flank10"].fillna("").astype(str).apply(
        lambda s: mask_telomere_fragments(s, motifs=TELO_MOTIFS, min_units=2)
    )
    out["downstream_flank10"] = out["downstream_flank10"].fillna("").astype(str).apply(
        lambda s: mask_telomere_fragments(s, motifs=TELO_MOTIFS, min_units=2)
    )
    return out


def build_logo_inputs(sub, unit_len):
    """Return (five_prime_strings, three_prime_strings) for rows whose
    first_unit AND last_unit both have length == unit_len.
    """
    sub = sub[(sub.first_unit != "") & (sub.last_unit != "")]
    sub = sub[sub.first_unit.str.len() == unit_len]
    sub = sub[sub.last_unit.str.len() == unit_len]
    if sub.empty:
        return [], []
    five = (sub.upstream_flank10 + sub.first_unit).tolist()
    three = (sub.last_unit + sub.downstream_flank10).tolist()
    return five, three


def plot_one(label, sub, outdir, background, flank_len=10, out_name=None,
             min_n=MIN_N_FOR_LOGO):
    """Build one figure per unit-length present in `sub` (multi-motif support).
    Returns the number of panels written.
    """
    if sub.empty:
        return 0
    # Per-unit-length stratification (so TTAGGG=6 + TTTAGGG=7 + TTAGG=5 each
    # get their own panel rather than the modal length silently dropping the
    # rest -- telotron_within_locus_motif_mixing_2026-06-03).
    unit_lens = sorted({len(u) for u in sub.first_unit if u}
                       | {len(u) for u in sub.last_unit if u})
    n_panels = 0
    for unit_len in unit_lens:
        five, three = build_logo_inputs(sub, unit_len)
        n = len(five)
        if n < min_n:
            print(f"  skip {label} unit={unit_len}: n={n} < {min_n}",
                  flush=True)
            continue
        fig, axes = plt.subplots(
            2, 1, figsize=(0.45 * (flank_len + unit_len) + 2.5, 4.6),
            layout="constrained",
        )
        fig.suptitle(
            f"{label}   n={n}   unit={unit_len} bp   "
            f"(flank masked; per-genome ACGT bg; e(n) corrected)",
            fontsize=10, fontweight="bold",
        )
        draw_panel(axes[0], five, "5p", flank_len, unit_len,
                   "5' end: upstream flank (left of |) + first repeat unit",
                   background=background)
        draw_panel(axes[1], three, "3p", flank_len, unit_len,
                   "3' end: last repeat unit + downstream flank (right of |)",
                   background=background)
        suffix = f"_unit{unit_len}"
        base = _slug(out_name or label) + suffix
        out = os.path.join(outdir, f"{base}.png")
        fig.savefig(out, dpi=160)
        plt.close(fig)
        n_panels += 1
    return n_panels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arrays", required=True,
                    help="work/results/interstitial_arrays.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank-len", type=int, default=10)
    ap.add_argument("--min-arrays", type=int, default=MIN_N_FOR_LOGO,
                    help="skip species with fewer arrays than this "
                         "(WebLogo small-sample threshold)")
    ap.add_argument("--filter-col", default=None)
    ap.add_argument("--filter-value", default=None)
    ap.add_argument("--refseq-dir", default="data/raw/refseq",
                    help="root for RefSeq FASTAs (for per-genome ACGT bg)")
    ap.add_argument("--tara-dir", default="data/raw/tara",
                    help="root for Tara MAG FASTAs (for per-genome ACGT bg)")
    ap.add_argument("--min-its-units", type=int, default=4,
                    help="Giulotto gold-standard minimum tandem-unit count")
    ap.add_argument("--min-its-ident", type=float, default=6.0 / 7.0,
                    help="Giulotto gold-standard minimum tandem identity")
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
        print(f"filter {args.filter_col} in {vals}: {before} -> {len(df)}",
              flush=True)

    # Gold-standard ITS filter (Giulotto): drop the ~73% of rows that fail
    # >=4 units / <1mm/unit. its_detection_gold_standard.
    df, (n_in, n_pass) = apply_gold_standard(
        df, min_units=args.min_its_units, min_ident=args.min_its_ident,
    )
    print(f"ITS gold standard (>={args.min_its_units} units, "
          f"ident>={args.min_its_ident:.3f}): {n_in} -> {n_pass}", flush=True)

    # Mask telomere rotations + revcomps + tandem >=2 runs out of flanks BEFORE
    # any logo construction. eimeria_subtelomere_rlfs_2026-05-15 retraction.
    df = mask_flanks(df, args.flank_len)

    # ----- overall (all species) — uniform background fallback -----
    plot_one("ALL_species", df, args.outdir, background=UNIFORM_BG,
             flank_len=args.flank_len, min_n=args.min_arrays)
    print(f"wrote ALL_species logos  (n={len(df)})", flush=True)

    # ----- per species: per-genome ACGT background -----
    bg_cache = {}
    for (gid, organism), sub in df.groupby(["genome_id", "organism"], sort=False):
        if len(sub) < args.min_arrays:
            continue
        # Per-genome ACGT background; falls back to uniform if FASTA missing.
        if gid not in bg_cache:
            try:
                fa = find_genome_fasta(
                    gid, args.refseq_dir, args.tara_dir, required=True,
                )
                bg_cache[gid] = compute_acgt_background(fa)
            except FileNotFoundError:
                print(f"  warn {gid}: no FASTA -> uniform background",
                      flush=True)
                bg_cache[gid] = dict(UNIFORM_BG)
        bg = bg_cache[gid]
        label = f"{organism} ({gid})"
        out_name = f"{_slug(gid)}_{_slug(organism)}"
        n_panels = plot_one(label, sub, args.outdir, background=bg,
                            flank_len=args.flank_len, out_name=out_name,
                            min_n=args.min_arrays)
        if n_panels:
            print(f"  {gid}: n={len(sub)}  panels={n_panels}  "
                  f"bg={ {k: round(v, 3) for k, v in bg.items()} }",
                  flush=True)


if __name__ == "__main__":
    main()
