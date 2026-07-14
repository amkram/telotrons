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
        try:
            for r in sub.itertuples(index=False):
                chrom = seqs.get(r.seqid, "")
                if not chrom:
                    continue
                L = len(chrom)
                # Requested flank window (may extend off contig); track actual+pad
                ls_req = r.start - 1 - args.flank
                re_req = r.end + args.flank
                ls = max(0, ls_req)
                le = r.start - 1
                rs = r.end
                re_ = min(L, re_req)
                left_pad_n = max(0, -ls_req)            # N's prepended to left flank (genomic)
                right_pad_n = max(0, re_req - L)        # N's appended to right flank (genomic)
                # Genomic-orientation flanks before strand flip + padding
                left_g = ("N" * left_pad_n) + chrom[ls:le]
                right_g = chrom[rs:re_] + ("N" * right_pad_n)
                # Build all sequences in display (spliced) orientation
                genomic = chrom[r.start - 1:r.end]
                if r.strand == "-":
                    intron = rc(genomic)
                    # On minus strand, display "left" is rc of genomic-right (and vice versa)
                    left = rc(right_g)
                    right = rc(left_g)
                    # Display-orientation pad counts (rc swaps + reverses)
                    disp_left_pad = right_pad_n
                    disp_right_pad = left_pad_n
                else:
                    intron = genomic
                    left = left_g
                    right = right_g
                    disp_left_pad = left_pad_n
                    disp_right_pad = right_pad_n

                arch = getattr(r, "architecture", "") if has_arch else ""
                linker_seq = getattr(r, "linker_seq", "") if has_arch else ""
                is_linker_arch = isinstance(linker_seq, str) and linker_seq and "linker" in arch

                if is_linker_arch:
                    # Prefer classifier-emitted linker_start/linker_end (0-based, into
                    # the spliced-orientation intron) — see classify_telotron_architecture.py.
                    # Falls back to string-search only if those columns are absent or unparseable.
                    lstart = getattr(r, "linker_start", None) if has_arch else None
                    lend = getattr(r, "linker_end", None) if has_arch else None
                    try:
                        lstart_i = int(lstart)
                        lend_i = int(lend)
                        if 0 <= lstart_i < lend_i <= len(intron):
                            arr1 = intron[:lstart_i]
                            lnk = intron[lstart_i:lend_i]
                            arr2 = intron[lend_i:]
                        else:
                            raise ValueError
                    except (TypeError, ValueError):
                        pos = intron.find(linker_seq)
                        if pos != -1:
                            arr1 = intron[:pos]
                            lnk = intron[pos:pos + len(linker_seq)]
                            arr2 = intron[pos + len(linker_seq):]
                        else:
                            arr1, lnk, arr2 = intron, "", ""
                else:
                    arr1, lnk, arr2 = intron, "", ""

                # Record actual emitted flank length AND N-pad counts so downstream
                # tools can detect contig-edge truncation (review finding 3).
                header = (f">{gid}|{r.seqid}|{r.start}-{r.end}|strand={r.strand}"
                          f"|tx={r.tx_id}|gene={r.gene_id}|len={r.intron_len}"
                          f"|motif={r.motif}|arch={arch}"
                          f"|left_len={len(left)}|right_len={len(right)}"
                          f"|left_pad_n={disp_left_pad}|right_pad_n={disp_right_pad}")
                intron_out = arr1 + lnk + arr2

                arch_key = arch if isinstance(arch, str) and arch else "Unknown"
                fa_out, fl_out = _get_handles(arch_key)
                fa_out.write(header + "\n" + intron_out + "\n")
                if is_linker_arch and lnk:
                    fl_out.write(header + "\n" + left + " " + arr1 + " " + lnk + " " + arr2 + " " + right + "\n")
                else:
                    fl_out.write(header + "\n" + left + " " + intron_out + " " + right + "\n")
                n += 1
        finally:
            for h in fa_handles.values():
                try:
                    h.close()
                except Exception:
                    pass
            for h in fl_handles.values():
                try:
                    h.close()
                except Exception:
                    pass
        per_arch = sub.groupby("architecture").size().to_dict() if "architecture" in sub.columns else {}
        print(f"  {gid}: {n} loci  " + " ".join(f"{a}={c}" for a, c in per_arch.items()))


if __name__ == "__main__":
    main()
