#!/usr/bin/env python
"""Rebuild the Eimeria subtelomere reference with Telogator2-compatible naming.

Telogator2 requires contig names of the form <sample>_chr<N><p|q> (single
underscore total), because telogator2.py:1051-1052 split on '_' to recover
the sample / chr tokens, and source/tg_util.py:LEXICO_2_IND keys are
'chr1'..'chr22'/'chrX'... etc. Our T1+T2 outputs use
<sample>_<accession>_<arm> (e.g. Eacervulina_OZ414021.1_p) which crashes
the final-allele sort.

This script renames the existing all_subtelomeres.fasta in-place and writes
a per-arm rename map.

Chromosome numbering: each Eimeria T2T has 15 nuclear contigs (the memory's
N=14 vs. N=15 puzzle: each species carries one ~1 Mb micro-chromosome that
is bi-capped and therefore a real micro-chromosome). We assign chr1..chr15
in descending order of contig size, which matches the published assembly
convention (largest=chr1).

LEXICO_2_IND in source/tg_util.py only ships chr1..chr22 + chrX/Y/M/F/B/U.
chr15 is present; we are safe up to chr22. Eimeria has 15 nuclear contigs.

Inputs:
  - cap_survey.tsv
  - subtelomeres/all_subtelomeres.fasta (the T2 deliverable)
Outputs (under work/results/telogator2_ref/):
  - all_subtelomeres.telogator2.fasta  (renamed FASTA, ready for `-t`)
  - rename_map.tsv                     (old_name <TAB> new_name)
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

REF_DIR = Path("work/results/telogator2_ref")
CAP_SURVEY = REF_DIR / "cap_survey.tsv"
IN_FASTA = REF_DIR / "subtelomeres" / "all_subtelomeres.fasta"
OUT_FASTA = REF_DIR / "all_subtelomeres.telogator2.fasta"
RENAME_MAP = REF_DIR / "rename_map.tsv"


def load_nuclear_by_species() -> dict[str, list[tuple[str, int]]]:
    by_species: dict[str, list[tuple[str, int]]] = {}
    with CAP_SURVEY.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("molecule_type") != "nuclear":
                continue
            by_species.setdefault(r["species"], []).append((r["contig_id"], int(r["length"])))
    return by_species


def build_rename(by_species: dict[str, list[tuple[str, int]]]) -> dict[str, str]:
    """Assign chr1..chrN by descending length per species; return
    {<sample>_<accession>_<p|q>: <sample>_chr<N><p|q>}."""
    rename: dict[str, str] = {}
    for species, contigs in by_species.items():
        contigs_sorted = sorted(contigs, key=lambda x: -x[1])
        for idx, (contig_id, _length) in enumerate(contigs_sorted, start=1):
            chrom = f"chr{idx}"
            for arm in ("p", "q"):
                old = f"{species}_{contig_id}_{arm}"
                new = f"{species}_{chrom}{arm}"
                rename[old] = new
    return rename


def parse_fasta(path: Path):
    name = None
    seq_lines: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_lines)
                name = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
    if name is not None:
        yield name, "".join(seq_lines)


def main() -> int:
    by_species = load_nuclear_by_species()
    print(f"Loaded cap-survey nuclear contigs: " +
          ", ".join(f"{sp}={len(v)}" for sp, v in by_species.items()))

    rename = build_rename(by_species)
    print(f"Built {len(rename)} potential arm -> chr renames "
          f"({len(rename)//2} arms x 2 arm-suffixes)")

    n_written = 0
    n_skipped = 0
    used_map: dict[str, str] = {}
    with OUT_FASTA.open("w") as out_fh:
        for name, seq in parse_fasta(IN_FASTA):
            new = rename.get(name)
            if new is None:
                print(f"  [warn] no rename for {name}; skipping", file=sys.stderr)
                n_skipped += 1
                continue
            used_map[name] = new
            out_fh.write(f">{new}\n")
            for i in range(0, len(seq), 80):
                out_fh.write(seq[i:i+80] + "\n")
            n_written += 1

    with RENAME_MAP.open("w") as fh:
        fh.write("old_name\tnew_name\n")
        for old, new in sorted(used_map.items()):
            fh.write(f"{old}\t{new}\n")

    print(f"Wrote {n_written} contigs to {OUT_FASTA} (skipped {n_skipped})")
    print(f"Wrote rename map -> {RENAME_MAP}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
