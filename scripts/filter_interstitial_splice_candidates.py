#!/usr/bin/env python3
"""Filter interstitial arrays to putative-telotron splice candidates.

Keep only arrays where:
  - the array is C-rich on the PLUS strand (i.e. `strand == "-"` in the
    interstitial output, which is the strand label we assign when the
    plus-strand sequence had to be reverse-complemented to reach the
    canonical G-rich telomeric orientation), AND
  - the 2 bp immediately upstream of the array on the PLUS strand is "GT"
    (a possible 5' donor splice site for a plus-strand intron), OR
  - the 2 bp immediately downstream of the array on the PLUS strand is "AG"
    (a possible 3' acceptor splice site).

These are loci that look like unannotated telotrons sitting inside a
plus-strand gene's intron, where the array reads CCCTAA on the mRNA strand.
"""
import argparse

import pandas as pd

from _common import rc


def plus_strand_2bp(row):
    """Return (upstream_2_plus, downstream_2_plus) on the genome plus strand.

    Stored flanks are in canonical (G-rich) telomeric orientation:
      strand == "+":  upstream_flank10 / downstream_flank10 are already on +
      strand == "-":  they're RC'd; swap and rc back to recover plus context.
    """
    up_canon  = str(row.upstream_flank10 or "")
    dn_canon  = str(row.downstream_flank10 or "")
    if row.strand == "+":
        upstream_2_plus   = up_canon[-2:] if len(up_canon)   >= 2 else ""
        downstream_2_plus = dn_canon[:2]  if len(dn_canon)   >= 2 else ""
    else:
        # canonical_upstream_flank10 = rc(chrom[end:end+10])  → plus downstream
        # canonical_downstream_flank10 = rc(chrom[start-10:start]) → plus upstream
        upstream_2_plus   = rc(dn_canon[:2])  if len(dn_canon)  >= 2 else ""
        downstream_2_plus = rc(up_canon[-2:]) if len(up_canon) >= 2 else ""
    return upstream_2_plus, downstream_2_plus


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arrays", required=True,
                    help="results/interstitial_arrays.tsv")
    ap.add_argument("--out", required=True,
                    help="filtered TSV (same schema + 3 derived columns)")
    args = ap.parse_args()

    df = pd.read_csv(args.arrays, sep="\t", low_memory=False)
    df = df.fillna({"upstream_flank10": "", "downstream_flank10": ""})

    up2, dn2 = [], []
    for r in df.itertuples(index=False):
        u, d = plus_strand_2bp(r)
        up2.append(u); dn2.append(d)
    df["plus_upstream_2bp"]   = up2
    df["plus_downstream_2bp"] = dn2

    keep = (
        (df.strand.astype(str) == "-")
        & ((df.plus_upstream_2bp.str.upper()   == "GT")
           | (df.plus_downstream_2bp.str.upper() == "AG"))
    )
    out = df[keep].copy()
    out["splice_motif_match"] = (
        out.plus_upstream_2bp.str.upper().eq("GT").map({True: "GT_donor", False: ""})
        + out.plus_downstream_2bp.str.upper().eq("AG").map({True: "_AG_acceptor", False: ""})
    ).str.strip("_")
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}  "
          f"({len(out)} / {len(df)} C-rich arrays with GT or AG splice signal "
          f"across {out.genome_id.nunique()} species)")


if __name__ == "__main__":
    main()
