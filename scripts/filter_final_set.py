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
    ap.add_argument("--include-partial", action="store_true",
                    help="Admit a flagged 'partial' tier (admission_pathway='partial'): a single-arm "
                         "clean array (>= --partial-min-hits tandem units on the dominant strand) at "
                         "moderate fraction (>= --partial-min-repeat-frac, below the strict single-array "
                         "threshold and not bidirectional). Captures partial / intermediate-state "
                         "telotrons the strict rules silently drop; flagged so downstream can stratify "
                         "them out of the strict counts.")
    ap.add_argument("--partial-min-repeat-frac", type=float, default=0.40)
    ap.add_argument("--partial-min-hits", type=int, default=3,
                    help="Min motif hits on the dominant strand for a partial array (excludes "
                         "scattered-hit noise).")
    ap.add_argument("--require-canonical-splice", action="store_true")
    ap.add_argument("--collapse-unique-loci", action="store_true")
    ap.add_argument("--require-terminal-motif-match", action="store_true",
                    help="Intron motif must equal the genome's dominant telomere motif. "
                         "CIRCULAR for downstream motif comparisons: the admitted set is by "
                         "construction motif-matched, so any figure that asks 'do telotron "
                         "motifs match host telomeres?' is definitional. When on, a "
                         "companion `*_no_motif_gate.tsv` is written as a sensitivity view.")
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
    # Partial tier (opt-in via --include-partial): a single-arm clean array
    # (>= partial_min_hits tandem units on one strand) at moderate fraction that meets
    # NEITHER the strict single-array (frac>=min_repeat_frac) NOR the bidirectional rule.
    # Captures partial/intermediate-state telotrons (a short array in a longer, partly-
    # ancestral intron) that the strict rules drop; the hit-count gate excludes the
    # scattered-telomeric-hit noise that dominates the 0.40-0.85 frac band.
    if args.include_partial:
        partial_pass = (
            ~single_pass & ~bidir_pass
            & (loci.telomeric_frac >= args.partial_min_repeat_frac)
            & ((loci.fwd_hits >= args.partial_min_hits) | (loci.rev_hits >= args.partial_min_hits))
        )
    else:
        partial_pass = pd.Series(False, index=loci.index)
    # Record which admission rule each locus passed, so downstream analyses can
    # stratify the strict single-array set (telomeric_frac >= min_repeat_frac), the
    # looser bidirectional set, and the flagged partial tier. Most loci enter bidir.
    loci = loci.assign(admission_pathway=np.where(
        single_pass & bidir_pass, "both",
        np.where(single_pass, "single_array",
        np.where(bidir_pass, "bidirectional",
        np.where(partial_pass, "partial", "none")))))
    admitted = loci[single_pass | bidir_pass | partial_pass].copy()

    def _apply_post_gate_filters(df):
        """canonical-splice filter + collapse_unique_loci, isolated so the
        `_no_motif_gate` companion can be produced with the same schema."""
        out = df
        if args.require_canonical_splice:
            out = out[out.splice_class == "GT-AG"]
        if args.collapse_unique_loci:
            # Collapse alt-spliced paralogs sharing (genome_id, seqid, start, end, strand)
            # into a single representative row, but PRESERVE tx_ids/gene_ids of the
            # dropped paralogs as comma-joined columns so downstream host-gene analyses
            # (telotron_ortholog_align, etc.) can recover alt-splice incidence.
            dup_keys = ["genome_id", "seqid", "start", "end", "strand"]
            if {"tx_id", "gene_id"}.issubset(out.columns):
                agg_ids = (out.groupby(dup_keys, dropna=False)
                              .agg(tx_ids=("tx_id", lambda s: ",".join(sorted(set(str(v) for v in s if pd.notna(v))))),
                                   gene_ids=("gene_id", lambda s: ",".join(sorted(set(str(v) for v in s if pd.notna(v))))))
                              .reset_index())
            else:
                agg_ids = None
            out = (out.sort_values("telomeric_frac", ascending=False)
                      .drop_duplicates(dup_keys))
            if agg_ids is not None:
                out = out.merge(agg_ids, on=dup_keys, how="left")
        return out

    # Sensitivity companion: apply the same post-gate filters WITHOUT the
    # (circular) motif gate. Same schema as final.tsv so downstream figures
    # can consume either without special-casing.
    final_no_motif_gate = _apply_post_gate_filters(admitted).copy()

    final = admitted
    if args.require_terminal_motif_match and len(final):
        final = final[final.terminal_motif.notna() & (final.terminal_motif != "")
                      & (final.motif == final.terminal_motif)]
    final = _apply_post_gate_filters(final)

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
    # Disambiguate true negatives from species whose candidates were killed by
    # require_terminal_motif_match / require_canonical_splice / collapse. The scan
    # summary's `telotron_candidates` is the pre-filter candidate count.
    # FIXES_PRIORITIZED.md unnumbered (review_survey_core.md MEDIUM, line 69).
    if "telotron_candidates" in species.columns:
        species["n_pre_filter_candidates"] = species["telotron_candidates"].fillna(0).astype(int)
        # Only species with zero pre-filter candidates AND a successful scan
        # are genuine negatives. Exclude rows tagged NO_FASTA / ERROR:* — those
        # never got scanned and mixing them into negatives inflates the
        # "true-negative" pool with missing-input artifacts.
        status = species.get("status", pd.Series([""] * len(species))).fillna("").astype(str)
        scanned_ok = ~status.str.startswith(("NO_FASTA", "NO_GFF", "ERROR"))
        negatives = species[(species.telotrons == 0) & (species.n_pre_filter_candidates == 0) & scanned_ok].copy()
    else:
        negatives = species[species.telotrons == 0].copy()

    final.to_csv(args.final, sep="\t", index=False)
    species.to_csv(args.species, sep="\t", index=False)
    negatives.to_csv(args.negatives, sep="\t", index=False)
    # Companion sensitivity output (only meaningful when --require-terminal-motif-match
    # is set): the pre-motif-gate table.
    if args.require_terminal_motif_match:
        final_no_motif_gate.to_csv(
            args.final.replace(".tsv", "_no_motif_gate.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
