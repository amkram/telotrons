#!/usr/bin/env python3
"""Lineage-specific substitution scan on a TERT ortholog alignment.

Goal: find positions that are conserved across telotron-POOR coccidia (and
broadly across distant outgroups) but specifically and consistently replaced
in the telotron-RICH Eimeria clade -- candidate functional changes that could
relate to telomerase promiscuity. Also reports per-domain Eimeria-vs-sister
divergence (relaxed-constraint proxy) using a reference protein's TRBD/RT
coordinates.

Sequence names must be prefixed EIM__ / POOR__ / OUT__.
"""
import argparse
from collections import Counter, OrderedDict


def read_aln(path):
    seqs = OrderedDict()
    h = None
    for l in open(path):
        l = l.rstrip()
        if l.startswith(">"):
            h = l[1:]
            seqs[h] = []
        else:
            seqs[h].append(l)
    return {k: "".join(v) for k, v in seqs.items()}


def consensus(col, names, minfrac):
    aa = [col[n] for n in names if col[n] != "-"]
    if not aa:
        return None, 0.0, 0
    c = Counter(aa)
    res, cnt = c.most_common(1)[0]
    return res, cnt / len(names), len(aa)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True)
    ap.add_argument("--ref", default="POOR__Cyclospora",
                    help="reference seq for ungapped coordinate mapping")
    ap.add_argument("--out", required=True)
    ap.add_argument("--eim-frac", type=float, default=0.8,
                    help="min fraction of Eimeria sharing the derived residue")
    ap.add_argument("--poor-frac", type=float, default=0.8,
                    help="min fraction of telotron-poor coccidia sharing the ancestral residue")
    args = ap.parse_args()

    aln = read_aln(args.aln)
    names = list(aln)
    L = len(next(iter(aln.values())))
    EIM = [n for n in names if n.startswith("EIM__")]
    POOR = [n for n in names if n.startswith("POOR__")]
    OUT = [n for n in names if n.startswith("OUT__")]

    # reference ungapped coordinate at each column
    ref = aln.get(args.ref) or aln[POOR[0]]
    refpos = []
    p = 0
    for ch in ref:
        if ch != "-":
            p += 1
        refpos.append(p)

    hits = []
    for i in range(L):
        col = {n: aln[n][i] for n in names}
        eim_res, eim_frac, eim_n = consensus(col, EIM, args.eim_frac)
        poor_res, poor_frac, poor_n = consensus(col, POOR, args.poor_frac)
        out_res, out_frac, out_n = consensus(col, OUT, 0.0)
        if eim_res is None or poor_res is None:
            continue
        # Eimeria-specific replacement: poor coccidia conserved (ancestral),
        # outgroups agree with poor (so ancestral state is poor's), Eimeria
        # consistently different.
        if (poor_frac >= args.poor_frac and eim_frac >= args.eim_frac
                and eim_res != poor_res
                and eim_n >= 4 and poor_n >= 4):
            ancestral_shared = (out_res == poor_res)  # outgroup confirms ancestral
            hits.append({
                "aln_col": i + 1, "ref_pos": refpos[i],
                "eimeria_res": eim_res, "eimeria_frac": round(eim_frac, 2),
                "poor_res": poor_res, "poor_frac": round(poor_frac, 2),
                "outgroup_res": out_res, "outgroup_frac": round(out_frac, 2),
                "outgroup_confirms_ancestral": ancestral_shared,
            })

    # per-residue identity of each Eimeria to the POOR consensus (divergence proxy)
    poor_cons = []
    for i in range(L):
        col = {n: aln[n][i] for n in POOR}
        r, f, n = consensus(col, POOR, 0.0)
        poor_cons.append(r if (r and f >= 0.6) else None)

    def pct_id_to_poor(seqname):
        s = aln[seqname]
        match = tot = 0
        for i in range(L):
            if poor_cons[i] is None or s[i] == "-":
                continue
            tot += 1
            if s[i] == poor_cons[i]:
                match += 1
        return 100 * match / tot if tot else 0.0

    with open(args.out, "w") as fh:
        fh.write("# Eimeria-specific substitutions at positions conserved in telotron-poor coccidia\n")
        fh.write(f"# Eimeria seqs: {len(EIM)}  telotron-poor coccidia: {len(POOR)}  outgroups: {len(OUT)}\n")
        fh.write(f"# reference (coords): {args.ref}\n")
        fh.write(f"# total Eimeria-specific replacements: {len(hits)} "
                 f"({sum(h['outgroup_confirms_ancestral'] for h in hits)} with outgroup-confirmed ancestral state)\n#\n")
        fh.write("aln_col\tref_pos\teimeria_res\teimeria_frac\tpoor_res\tpoor_frac\toutgroup_res\toutgroup_confirms_ancestral\n")
        for h in hits:
            fh.write(f"{h['aln_col']}\t{h['ref_pos']}\t{h['eimeria_res']}\t{h['eimeria_frac']}\t"
                     f"{h['poor_res']}\t{h['poor_frac']}\t{h['outgroup_res']}\t{h['outgroup_confirms_ancestral']}\n")

    print(f"Eimeria-specific replacements (conserved-in-sisters -> changed-in-Eimeria): {len(hits)}")
    print(f"  with outgroup-confirmed ancestral state (high confidence): "
          f"{sum(h['outgroup_confirms_ancestral'] for h in hits)}")
    print("\n%identity of each sequence to the telotron-poor-coccidia consensus:")
    for n in EIM + POOR:
        print(f"  {n:24s} {pct_id_to_poor(n):5.1f}%")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
