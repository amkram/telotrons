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


def stack(seqs, length):
    counts = {b: [0] * length for b in ACGTN}
    for s in seqs:
        if len(s) != length:
            continue
        for i, ch in enumerate(s.upper()):
            if ch in counts:
                counts[ch][i] += 1
    df = pd.DataFrame(counts)
    row_sums = df.sum(axis=1).replace(0, 1)
    return df.div(row_sums, axis=0)


def to_info(freq):
    return logomaker.transform_matrix(freq, from_type="probability",
                                      to_type="information")


def draw(ax, seqs, length, title, signal_pos):
    info = to_info(stack(seqs, length))
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


def plot_one(label, donor_seqs, acceptor_seqs, outdir):
    if not donor_seqs and not acceptor_seqs:
        return False
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 4.6), layout="constrained")
    n = max(len(donor_seqs), len(acceptor_seqs))
    fig.suptitle(f"{label}   n={n}   ±10 bp of splice signals",
                 fontsize=11, fontweight="bold")
    # donor: positions 0..9 = exon, 10..11 = GT, 12..21 = intron
    draw(axes[0], donor_seqs, 22, "5' splice site (donor)  exon | GT | intron",
         signal_pos=10)
    # acceptor: 0..9 = intron, 10..11 = AG, 12..21 = exon
    draw(axes[1], acceptor_seqs, 22, "3' splice site (acceptor)  intron | AG | exon",
         signal_pos=10)
    out = os.path.join(outdir, f"{_slug(label)}.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flanked-dir", required=True,
                    help="root dir of per-species *.txt flanked files")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-records", type=int, default=5,
                    help="skip species with fewer than this many records")
    ap.add_argument("--annotation-tsv", default=None,
                    help="TSV with genome_id/seqid/start/end + filter column")
    ap.add_argument("--filter-col", default=None)
    ap.add_argument("--filter-value", default=None)
    args = ap.parse_args()

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

    plot_one("ALL_species", all_donor, all_acceptor, args.outdir)
    print(f"wrote ALL_species  (n={len(all_donor)})", flush=True)
    for sp, (d, a) in sorted(by_species.items()):
        if len(d) < args.min_records:
            continue
        plot_one(sp, d, a, args.outdir)
        print(f"  {sp}: n={len(d)}", flush=True)


if __name__ == "__main__":
    main()
