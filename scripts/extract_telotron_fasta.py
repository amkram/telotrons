#!/usr/bin/env python3
"""Per-species, per-architecture FASTA + flanked-text extracts of filtered telotron loci.

Outputs (one subdir per positive species, one file per architecture):
    <fasta_dir>/<gid>_<organism>/<arch>.fa
    <flanked_dir>/<gid>_<organism>/<arch>.txt

Flanked sequence line:
    non-linker arch:  [LEFT100] [INTRON] [RIGHT100]
    linker arch:      [LEFT100] [ARRAY1] [LINKER] [ARRAY2] [RIGHT100]

Architectures (see classify_telotron_architecture.py for the assignment rules):
    GT-F-AG, GT-R-AG, GT-F-R-AG, GT-F-linker-R-AG, GT-R-linker-F-AG, Other, Unknown
"""
import argparse
import os

import pandas as pd

from _common import rc, slug as _slug, load_fasta, find_genome_fasta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True, help="final_telotron_set_architecture.tsv")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--flanked-dir", required=True)
    ap.add_argument("--flank", type=int, default=100)
    args = ap.parse_args()

    os.makedirs(args.fasta_dir, exist_ok=True)
    os.makedirs(args.flanked_dir, exist_ok=True)

    df = pd.read_csv(args.final, sep="\t")
    has_arch = "architecture" in df.columns

    for gid, sub in df.groupby("genome_id", sort=False):
        try:
            fa = find_genome_fasta(gid, args.refseq_dir, args.tara_dir)
        except FileNotFoundError:
            print(f"skip {gid}: FASTA not found")
            continue
        seqs = load_fasta(fa)
        organism = sub.organism.iloc[0]
        slug = f"{_slug(gid)}_{_slug(organism)}"
        fa_dir = f"{args.fasta_dir}/{slug}"
        fl_dir = f"{args.flanked_dir}/{slug}"
        os.makedirs(fa_dir, exist_ok=True)
        os.makedirs(fl_dir, exist_ok=True)

        fa_handles = {}
        fl_handles = {}

        def _get_handles(arch_key):
            if arch_key not in fa_handles:
                fa_handles[arch_key] = open(f"{fa_dir}/{arch_key}.fa", "w")
                fl_handles[arch_key] = open(f"{fl_dir}/{arch_key}.txt", "w")
            return fa_handles[arch_key], fl_handles[arch_key]

        n = 0
        for r in sub.itertuples(index=False):
            chrom = seqs.get(r.seqid, "")
            if not chrom:
                continue
            L = len(chrom)
            ls = max(0, r.start - 1 - args.flank)
            le = r.start - 1
            rs = r.end
            re_ = min(L, r.end + args.flank)
            # Build all sequences in display (spliced) orientation
            genomic = chrom[r.start - 1:r.end]
            if r.strand == "-":
                intron = rc(genomic)
                left = rc(chrom[rs:re_])
                right = rc(chrom[ls:le])
            else:
                intron = genomic
                left = chrom[ls:le]
                right = chrom[rs:re_]

            arch = getattr(r, "architecture", "") if has_arch else ""
            linker_seq = getattr(r, "linker_seq", "") if has_arch else ""
            is_linker_arch = isinstance(linker_seq, str) and linker_seq and "linker" in arch

            if is_linker_arch:
                # linker_seq was extracted from spliced (display) orientation by classify script
                pos = intron.find(linker_seq)
                if pos != -1:
                    arr1 = intron[:pos]
                    lnk = intron[pos:pos + len(linker_seq)]
                    arr2 = intron[pos + len(linker_seq):]
                else:
                    arr1, lnk, arr2 = intron, "", ""
            else:
                arr1, lnk, arr2 = intron, "", ""

            header = (f">{gid}|{r.seqid}|{r.start}-{r.end}|strand={r.strand}"
                      f"|tx={r.tx_id}|gene={r.gene_id}|len={r.intron_len}"
                      f"|motif={r.motif}|arch={arch}")
            intron_out = arr1 + lnk + arr2

            arch_key = arch if isinstance(arch, str) and arch else "Unknown"
            fa_out, fl_out = _get_handles(arch_key)
            fa_out.write(header + "\n" + intron_out + "\n")
            if is_linker_arch and lnk:
                fl_out.write(header + "\n" + left + " " + arr1 + " " + lnk + " " + arr2 + " " + right + "\n")
            else:
                fl_out.write(header + "\n" + left + " " + intron_out + " " + right + "\n")
            n += 1
        for h in fa_handles.values(): h.close()
        for h in fl_handles.values(): h.close()
        per_arch = sub.groupby("architecture").size().to_dict() if "architecture" in sub.columns else {}
        print(f"  {gid}: {n} loci  " + " ".join(f"{a}={c}" for a, c in per_arch.items()))


if __name__ == "__main__":
    main()
