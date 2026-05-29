#!/usr/bin/env python3
"""Filter candidate loci to the final telotron set; build species-level summary and negative controls."""
import argparse

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True)
    ap.add_argument("--introns", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--min-repeat-frac", type=float, default=0.85,
                    help="Single-array threshold: dense G- or C-rich introns must meet this.")
    ap.add_argument("--bidir-min-repeat-frac", type=float, default=0.40,
                    help="Looser threshold for bidirectional (diverging/converging) introns "
                         "that satisfy --bidir-min-hits on both strands.")
    ap.add_argument("--bidir-min-hits", type=int, default=3,
                    help="Per-strand hit count required to qualify as bidirectional.")
    ap.add_argument("--require-canonical-splice", action="store_true")
    ap.add_argument("--collapse-unique-loci", action="store_true")
    ap.add_argument("--require-terminal-motif-match", action="store_true",
                    help="Intron motif must equal the genome's dominant telomere motif.")
    ap.add_argument("--final", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--negatives", required=True)
    args = ap.parse_args()

    loci = pd.read_csv(args.loci, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")

    # Single-array: meets the strict fraction.
    # Bidirectional: looser fraction AND ≥N hits on EACH strand (diverging/converging signature).
    single_pass = loci.telomeric_frac >= args.min_repeat_frac
    bidir_pass = (
        (loci.telomeric_frac >= args.bidir_min_repeat_frac)
        & (loci.fwd_hits >= args.bidir_min_hits)
        & (loci.rev_hits >= args.bidir_min_hits)
    )
    # Record which admission rule each locus passed, so downstream analyses can
    # stratify the strict single-array set (telomeric_frac >= min_repeat_frac)
    # from the looser bidirectional set (frac >= bidir_min_repeat_frac with
    # >= bidir_min_hits on both strands). Most loci enter via the looser rule.
    loci = loci.assign(admission_pathway=np.where(
        single_pass & bidir_pass, "both",
        np.where(single_pass, "single_array", "bidirectional")))
    final = loci[single_pass | bidir_pass].copy()
    if args.require_terminal_motif_match and len(final):
        final = final[final.terminal_motif.notna() & (final.terminal_motif != "")
                      & (final.motif == final.terminal_motif)]
    if args.require_canonical_splice:
        final = final[final.splice_class == "GT-AG"]
    if args.collapse_unique_loci:
        final = (final.sort_values("telomeric_frac", ascending=False)
                      .drop_duplicates(["genome_id", "seqid", "start", "end", "strand"]))

    per_species = (final.groupby(["genome_id", "organism", "group", "source", "motif"], dropna=False)
                        .agg(telotrons=("seqid", "size"),
                             canonical_gt_ag=("splice_class", lambda x: (x == "GT-AG").sum()),
                             fwd_g_rich=("orientation", lambda x: (x == "Fwd/G-rich").sum()),
                             rev_c_rich=("orientation", lambda x: (x == "Rev/C-rich").sum()),
                             mixed=("orientation", lambda x: (x == "Mixed").sum()),
                             median_len=("intron_len", "median"))
                        .reset_index())

    species = summary.merge(per_species, on=["genome_id", "organism", "group", "source"], how="left")
    species["telotrons"] = species["telotrons"].fillna(0).astype(int)
    negatives = species[species.telotrons == 0].copy()

    final.to_csv(args.final, sep="\t", index=False)
    species.to_csv(args.species, sep="\t", index=False)
    negatives.to_csv(args.negatives, sep="\t", index=False)


if __name__ == "__main__":
    main()
