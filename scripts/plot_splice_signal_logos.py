#!/usr/bin/env python3
"""Sequence logos centered on each splice signal (donor GT and acceptor AG).

Input: a directory tree of `.txt` flanked files produced by
extract_telotron_fasta.py — each non-header line is space-separated:

    [LEFT100] [INTRON] [RIGHT100]          (non-linker arch)
    [LEFT100] [ARRAY1] [LINKER] [ARRAY2] [RIGHT100]  (linker arch)

For every record we slice two 22-bp windows in display (spliced) orientation:

    donor    = LEFT100[-10:]  + INTRON[:12]    (10 bp exon + GT + 10 bp intron)
    acceptor = INTRON[-12:]   + RIGHT100[:10]  (10 bp intron + AG + 10 bp exon)

For each species (and a combined "ALL_species") we draw a 2-row info-bit logo:
donor on top, acceptor on bottom. A dashed line marks the exon/intron boundary.
"""
import argparse
import glob
import os

import matplotlib.pyplot as plt
import pandas as pd
import logomaker

from _common import slug as _slug

ACGTN = list("ACGT")

# Schneider 1986 small-sample correction is applied automatically by
# logomaker only when from_type='counts'. Refuse to draw logos below this
# many records (verify_figures_visualization Finding 5; n>=20 keeps the
# e(n) bias term <~0.05 bits).
MIN_LOGO_N = 20


def load_background_freqs(path):
    """Load per-genome ACGT background frequencies.

    Expected TSV columns: genome_id, A, C, G, T (counts or probabilities).
    Returns dict[genome_id] -> {A,C,G,T: float (probabilities)}.
    """
    if not path or not os.path.exists(path):
        return {}
    bg = pd.read_csv(path, sep="\t")
    out = {}
    for _, r in bg.iterrows():
        try:
            vals = {b: float(r[b]) for b in ACGTN}
        except (KeyError, ValueError):
            continue
        total = sum(vals.values())
        if total <= 0:
            continue
        out[str(r["genome_id"])] = {b: vals[b] / total for b in ACGTN}
    return out


def background_for_species(sp_slug, by_genome):
    """Look up an ACGT background for a species slug.

    The flanked-dir uses species slugs; if the background table is keyed
    by genome_id and `sp_slug` matches an entry, return it; otherwise
    return None (caller falls back to uniform with a printed warning).
    """
    if not by_genome:
        return None
    if sp_slug in by_genome:
        return by_genome[sp_slug]
    # SPECIES_BY_ACC maps GCF_xxx -> species name; try the reverse via slug.
    try:
        from _common import SPECIES_BY_ACC, slug as _s
        for acc, species in SPECIES_BY_ACC.items():
            if _s(species) == sp_slug and acc in by_genome:
                return by_genome[acc]
    except Exception:
        pass
    return None


def parse_header(line):
    """Parse the >gid|seqid|start-end|... headers written by extract_telotron_fasta."""
    parts = line[1:].split("|")
    if len(parts) < 3:
        return None
    gid, seqid, span = parts[0], parts[1], parts[2]
    if "-" not in span:
        return None
    try:
        start, end = span.split("-")
        start, end = int(start), int(end)
    except ValueError:
        return None
    return (gid, seqid, start, end)


def parse_flanked_dir(root, allowed=None):
    """Yield (species_slug, donor_window, acceptor_window) for every record.

    If `allowed` is a set of (gid, seqid, start, end) tuples, skip records
    whose header isn't in the set.
    """
    for sp_dir in sorted(os.listdir(root)):
        sp_path = os.path.join(root, sp_dir)
        if not os.path.isdir(sp_path):
            continue
        for fpath in glob.glob(os.path.join(sp_path, "*.txt")):
            with open(fpath) as fh:
                header = None
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith(">"):
                        header = line
                        continue
                    if allowed is not None:
                        key = parse_header(header) if header else None
                        if key is None or key not in allowed:
                            continue
                    parts = line.split(" ")
                    if len(parts) < 3:
                        continue
                    left, right = parts[0], parts[-1]
                    # Middle parts concatenate to the spliced-orientation intron
                    # (1 part = non-linker; 3 parts = arr1 + linker + arr2).
                    middle = "".join(parts[1:-1])
                    if len(left) < 10 or len(right) < 10 or len(middle) < 12:
                        continue
                    donor    = left[-10:] + middle[:12]
                    acceptor = middle[-12:] + right[:10]
                    yield sp_dir, donor, acceptor


def stack_counts(seqs, length):
    """Return per-position ACGT counts DataFrame (integer-valued)."""
    counts = {b: [0] * length for b in ACGTN}
    for s in seqs:
        if len(s) != length:
            continue
        for i, ch in enumerate(s.upper()):
            if ch in counts:
                counts[ch][i] += 1
    return pd.DataFrame(counts)


def to_info(counts_df, background=None):
    """counts -> information using Schneider 1986 small-sample correction.

    `from_type='counts'` invokes logomaker's pseudocount path (1.0 by
    default) which provides the e(n) small-sample bias term. If
    `background` is a dict {A,C,G,T: prob}, it is converted to an array
    ordered to match `counts_df.columns` before being passed to logomaker
    (logomaker accepts array-like / DataFrame, not dict).
    """
    bg_arg = None
    if background is not None:
        bg_arg = [float(background[c]) for c in counts_df.columns]
    return logomaker.transform_matrix(
        counts_df,
        from_type="counts",
        to_type="information",
        background=bg_arg,
        pseudocount=1,
    )


def draw(ax, seqs, length, title, signal_pos, background=None):
    counts = stack_counts(seqs, length)
    info = to_info(counts, background=background)
    logomaker.Logo(info, ax=ax, color_scheme="classic",
                   show_spines=False, fade_below=0.0)
    ax.set_ylim(0, 2)
    ax.set_yticks([0, 1, 2])
    ax.set_ylabel("bits")
    ax.set_title(title, fontsize=10)
    # Mark the splice-signal dinucleotide (positions signal_pos, signal_pos+1)
    ax.axvline(signal_pos - 0.5, color="#888", lw=0.8, ls="--")
    ax.axvline(signal_pos + 1.5, color="#888", lw=0.8, ls="--")
    ax.set_xticks(range(length))
    # Label positions relative to the splice signal (-10..-1, GT, +1..+10)
    labels = []
    for i in range(length):
        d = i - signal_pos
        if d == 0 or d == 1:
            labels.append("")
        elif d < 0:
            labels.append(str(d))
        else:
            labels.append(f"+{d - 1}")
    ax.set_xticklabels(labels, fontsize=7)


def plot_one(label, donor_seqs, acceptor_seqs, outdir, background=None,
             min_n=MIN_LOGO_N):
    if not donor_seqs and not acceptor_seqs:
        return False
    n = max(len(donor_seqs), len(acceptor_seqs))
    if n < min_n:
        print(f"skip {label}: n={n} < min_n={min_n} "
              "(Schneider 1986 small-sample regime)", flush=True)
        return False
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 4.6), layout="constrained")
    bg_note = " bg=uniform" if background is None else " bg=genome"
    fig.suptitle(f"{label}   n={n}   pm10 bp of splice signals{bg_note}",
                 fontsize=11, fontweight="bold")
    # donor: positions 0..9 = exon, 10..11 = GT, 12..21 = intron
    draw(axes[0], donor_seqs, 22, "5' splice site (donor)  exon | GT | intron",
         signal_pos=10, background=background)
    # acceptor: 0..9 = intron, 10..11 = AG, 12..21 = exon
    draw(axes[1], acceptor_seqs, 22, "3' splice site (acceptor)  intron | AG | exon",
         signal_pos=10, background=background)
    out = os.path.join(outdir, f"{_slug(label)}.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flanked-dir", required=True,
                    help="root dir of per-species *.txt flanked files")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-records", type=int, default=MIN_LOGO_N,
                    help=("skip species with fewer than this many records "
                          "(default %d = Schneider 1986 small-sample floor)"
                          % MIN_LOGO_N))
    ap.add_argument("--annotation-tsv", default=None,
                    help="TSV with genome_id/seqid/start/end + filter column")
    ap.add_argument("--filter-col", default=None)
    ap.add_argument("--filter-value", default=None)
    ap.add_argument("--background-tsv", default=None,
                    help=("per-genome ACGT background TSV "
                          "(columns: genome_id, A, C, G, T) for Schneider "
                          "& Stephens 1990 information content. If absent, "
                          "uniform background is used and noted on the plot."))
    args = ap.parse_args()

    by_genome_bg = load_background_freqs(args.background_tsv)
    if args.background_tsv and not by_genome_bg:
        print(f"warning: background TSV {args.background_tsv} "
              "could not be loaded; falling back to uniform.", flush=True)

    allowed = None
    if args.annotation_tsv and args.filter_col and args.filter_value is not None:
        ann = pd.read_csv(args.annotation_tsv, sep="\t", low_memory=False)
        before = len(ann)
        vals = [v for v in str(args.filter_value).split(",") if v]
        ann = ann[ann[args.filter_col].astype(str).isin(vals)]
        allowed = set(zip(ann.genome_id.astype(str), ann.seqid.astype(str),
                          ann.start.astype(int), ann.end.astype(int)))
        print(f"filter {args.filter_col} in {vals}: "
              f"{before} → {len(ann)} rows, {len(allowed)} unique loci",
              flush=True)

    os.makedirs(args.outdir, exist_ok=True)
    by_species = {}
    all_donor, all_acceptor = [], []
    for sp, donor, acceptor in parse_flanked_dir(args.flanked_dir, allowed=allowed):
        by_species.setdefault(sp, ([], []))
        by_species[sp][0].append(donor)
        by_species[sp][1].append(acceptor)
        all_donor.append(donor)
        all_acceptor.append(acceptor)

    # Pooled "ALL_species" panel uses average background across genomes seen,
    # weighted by per-species record count. Falls back to uniform if no bg
    # table was supplied.
    pooled_bg = None
    if by_genome_bg:
        agg = {b: 0.0 for b in ACGTN}
        wsum = 0.0
        for sp, (d, _a) in by_species.items():
            bg = background_for_species(sp, by_genome_bg)
            if bg is None:
                continue
            w = float(len(d))
            for b in ACGTN:
                agg[b] += w * bg[b]
            wsum += w
        if wsum > 0:
            pooled_bg = {b: agg[b] / wsum for b in ACGTN}

    plot_one("ALL_species", all_donor, all_acceptor, args.outdir,
             background=pooled_bg, min_n=args.min_records)
    print(f"wrote ALL_species  (n={len(all_donor)}) "
          f"bg={'genome-weighted' if pooled_bg else 'uniform'}", flush=True)
    for sp, (d, a) in sorted(by_species.items()):
        if len(d) < args.min_records:
            print(f"  {sp}: n={len(d)} < min_records={args.min_records} - skipped",
                  flush=True)
            continue
        bg = background_for_species(sp, by_genome_bg)
        plot_one(sp, d, a, args.outdir, background=bg,
                 min_n=args.min_records)
        print(f"  {sp}: n={len(d)} bg={'genome' if bg else 'uniform'}",
              flush=True)


if __name__ == "__main__":
    main()
