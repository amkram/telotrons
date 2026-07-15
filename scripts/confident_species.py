#!/usr/bin/env python3
"""Emit a species-level confident-telotron-bearing set from filter_final output.

A species is admitted as a confident bearer when it either
  (a) has >= --min-n telotrons that pass filter_final, OR
  (b) has >= 1 bidirectional telotron (GT-F-R-AG, GT-R-linker-F-AG, GT-F-linker-R-AG),
      a distinctive telomerase-mediated architecture unlikely by chance.

Output TSV: genome_id, organism, source, n_telotrons, n_bidirectional, admission.
Downstream analyses (gene-class, expression, nucleosome, ortholog panels) key
off this file so new bearer species flow through automatically."""
import argparse
import csv
from collections import Counter

BIDIR_ARCHS = {"GT-F-R-AG", "GT-R-linker-F-AG", "GT-F-linker-R-AG"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", required=True, help="final_telotron_set_architecture.tsv")
    ap.add_argument("--manifest", required=True, help="work/manifests/all_genomes.tsv")
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    src = {r["genome_id"]: r.get("source", "")
           for r in csv.DictReader(open(args.manifest), delimiter="\t")}

    n = Counter()
    n_bidir = Counter()
    organism = {}
    for r in csv.DictReader(open(args.arch), delimiter="\t"):
        gid = r["genome_id"]
        n[gid] += 1
        if r.get("architecture", "") in BIDIR_ARCHS:
            n_bidir[gid] += 1
        organism.setdefault(gid, r.get("organism", ""))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "organism", "source", "n_telotrons",
                    "n_bidirectional", "admission"])
        for gid in sorted(n):
            passes = []
            if n[gid] >= args.min_n:
                passes.append(f"n>={args.min_n}")
            if n_bidir[gid] >= 1:
                passes.append("bidirectional")
            if not passes:
                continue
            w.writerow([gid, organism[gid], src.get(gid, ""),
                        n[gid], n_bidir[gid], ",".join(passes)])
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
