#!/usr/bin/env python3
"""Build a control set of non-telotron introns (telomeric_frac < threshold).

For each genome with at least one telotron, sample up to --n-per-species
introns whose `telomeric_frac` is below `--max-frac` (default 0.10). The
output TSV mirrors the schema of final_telotron_set_architecture.tsv so the
existing extract_telotron_fasta.py can ingest it: an `architecture` column is
added with the constant value "control", and `linker_seq` is empty.
"""
import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--introns", required=True, help="results/all_introns_scanned.tsv")
    ap.add_argument("--final", required=True,
                    help="results/final_telotron_set.tsv (defines positive species)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-frac", type=float, default=0.10,
                    help="upper bound on telomeric_frac for control introns")
    ap.add_argument("--n-per-species", type=int, default=5000,
                    help="random sample size per species (cap)")
    ap.add_argument("--seed", type=int, default=0xC047F0)
    args = ap.parse_args()

    pos_species = set(pd.read_csv(args.final, sep="\t", usecols=["genome_id"]).genome_id)
    introns = pd.read_csv(args.introns, sep="\t", low_memory=False)

    sub = introns[introns.genome_id.isin(pos_species)
                  & (introns.telomeric_frac.astype(float) < args.max_frac)]
    pieces = [g.sample(min(len(g), args.n_per_species), random_state=args.seed)
              for _, g in sub.groupby("genome_id", sort=False)]
    sampled = pd.concat(pieces, ignore_index=True) if pieces else sub.head(0).copy()
    sampled["architecture"] = "control"
    if "linker_seq" not in sampled.columns:
        sampled["linker_seq"] = ""
    sampled.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}  "
          f"({len(sampled)} controls across {sampled.genome_id.nunique()} species; "
          f"telomeric_frac < {args.max_frac})")


if __name__ == "__main__":
    main()
