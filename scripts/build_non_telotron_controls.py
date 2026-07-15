#!/usr/bin/env python3
"""Build a control set of non-telotron introns (telomeric_frac < threshold).

For each genome (by default only genomes with at least one telotron; pass
``--all-genomes`` to draw across every surveyed genome), sample up to
``--n-per-species`` introns whose ``telomeric_frac`` is below ``--max-frac``
(default 0.01). The output TSV mirrors the schema of
``final_telotron_set_architecture.tsv`` so the existing ``extract_telotron_fasta.py``
can ingest it: an ``architecture`` column is added with the constant value
``"control"`` and ``linker_seq`` is empty.

Notes
-----
- The positive (telotron) set is *not* GT-AG-filtered by default
  (``config.yaml:require_canonical_splice: false``). To avoid the
  splice-class mismatch flagged in review_sequence_extraction_and_msa
  (Finding 16, FIXES_PRIORITIZED #40), the control set also defaults
  to NOT filtering by splice class. Pass ``--require-gt-ag`` to opt
  back into the GT-AG-only behaviour (deprecated alias:
  ``--allow-noncanonical`` is now the default and is a no-op kept for
  backwards compatibility).
- ``max_frac`` default tightened 0.10 -> 0.01 to keep partial-telotron
  intermediates (0.1 <= tf < 0.45, see ``telotron_intermediate_states_2026-06-03``)
  out of the negative-control set.
- Per-group ``random_state`` is now derived from ``hash((seed, gid))`` so
  re-runs that add a new species do not reshuffle the sample drawn for
  existing species (FIXES_PRIORITIZED #58).
"""
import argparse

import hashlib

import pandas as pd


def _group_seed(base_seed: int, gid: str) -> int:
    """Stable per-group seed. SHA1-based (not Python's built-in hash()) so it
    doesn't drift across interpreter runs — PYTHONHASHSEED randomization would
    otherwise reshuffle the per-species control sample on every rerun."""
    return int(hashlib.sha1(f"{base_seed}:{gid}".encode()).hexdigest()[:8], 16)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--introns", required=True, help="work/results/all_introns_scanned.tsv")
    ap.add_argument("--final", required=True,
                    help="work/results/final_telotron_set.tsv (defines positive species)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-frac", type=float, default=0.01,
                    help="upper bound on telomeric_frac for control introns "
                         "(tightened from 0.10 to 0.01 to exclude partial-telotron "
                         "intermediates; see telotron_intermediate_states_2026-06-03)")
    ap.add_argument("--n-per-species", type=int, default=500,
                    help="random sample size per species (cap). 5000 inflated "
                         "control variance in downstream Fisher/Mann-Whitney tests; "
                         "default lowered to 500 (~10x positive max).")
    ap.add_argument("--seed", type=int, default=0xC047F0)
    ap.add_argument("--all-genomes", action="store_true",
                    help="sample controls across ALL surveyed genomes, not only "
                         "positive-telotron species. Use for cross-species "
                         "comparisons (avoids the circular framing flagged in "
                         "feedback_intronic_ITS_circular.md). Default keeps the "
                         "within-positive-species behaviour for matched-host tests.")
    ap.add_argument("--require-gt-ag", action="store_true",
                    help="keep only GT-AG control introns. Default: no splice-class "
                         "filter (matches the positive set; see Snakemake config "
                         "`require_canonical_splice: false`). Without this flag, "
                         "controls share the splice-class distribution of positives "
                         "so boundary-k-mer comparisons are not confounded by a "
                         "splice-class filter applied to only one arm.")
    # Backwards-compatible no-op: previously, --allow-noncanonical opted *out*
    # of a GT-AG-only default. The default is now "allow noncanonical", so
    # this flag is a no-op kept so existing Snakemake rules keep working.
    ap.add_argument("--allow-noncanonical", action="store_true",
                    help="DEPRECATED: now the default. Use --require-gt-ag to opt "
                         "back into GT-AG-only sampling.")
    args = ap.parse_args()

    pos_species = set(pd.read_csv(args.final, sep="\t", usecols=["genome_id"]).genome_id)
    introns = pd.read_csv(args.introns, sep="\t", low_memory=False)

    if args.all_genomes:
        species_mask = pd.Series(True, index=introns.index)
    else:
        species_mask = introns.genome_id.isin(pos_species)

    sub = introns[species_mask
                  & (introns.telomeric_frac.astype(float) < args.max_frac)]
    if args.require_gt_ag and "splice_class" in sub.columns:
        sub = sub[sub.splice_class == "GT-AG"]
    pieces = [g.sample(min(len(g), args.n_per_species),
                       random_state=_group_seed(args.seed, gid))
              for gid, g in sub.groupby("genome_id", sort=False)]
    sampled = pd.concat(pieces, ignore_index=True) if pieces else sub.head(0).copy()
    sampled["architecture"] = "control"
    if "linker_seq" not in sampled.columns:
        sampled["linker_seq"] = ""
    sampled.to_csv(args.out, sep="\t", index=False)
    scope = "all surveyed genomes" if args.all_genomes else "positive-telotron species"
    splice_note = "GT-AG only" if args.require_gt_ag else "no splice-class filter"
    print(f"wrote {args.out}  "
          f"({len(sampled)} controls across {sampled.genome_id.nunique()} species; "
          f"telomeric_frac < {args.max_frac}; scope={scope}; {splice_note})")


if __name__ == "__main__":
    main()
