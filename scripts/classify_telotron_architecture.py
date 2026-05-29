#!/usr/bin/env python3
"""Classify each telotron locus by F/R-array architecture and emit boundary 6-mers.

Architectures (5 explicit + Other):
    GT-F-AG            : pure forward (G-rich) array
    GT-R-AG            : pure reverse (C-rich) array
    GT-F-R-AG          : convergent — F then R, no real linker between
    GT-F-linker-R-AG   : F then non-telomeric linker then R
    GT-R-linker-F-AG   : R then non-telomeric linker then F
    Other              : anything else

A "real linker" must be:
    - length >= MIN_LINKER (15 bp default)
    - <30% bases covered by F/R motif variants allowing 1 mismatch
This kills degenerate-telomeric "linkers" (the trap that bit us previously).

Output:
    --out-loci : per-locus TSV with architecture, donor_6mer, acceptor_6mer, linker_seq, linker_len
    --out-kmers: long-form TSV (genome_id, architecture, side, kmer, count)
"""
import argparse

import pandas as pd

from _common import rc, rotations, load_fasta, find_genome_fasta

MIN_ARRAY = 18           # minimum length of an F or R array to count
MAX_ARRAY_GAP = 2        # gap (bp) below which two motif hits merge into one array
MIN_LINKER = 15          # minimum gap length to call a linker
MAX_LINKER_TELO = 0.30   # max fraction of linker bases inside any motif (1-mm)


def find_intervals(seq, variants):
    intervals = []
    for v in variants:
        k = len(v)
        i = seq.find(v)
        while i != -1:
            intervals.append((i, i + k))
            i = seq.find(v, i + 1)
    intervals.sort()
    return intervals


def merge_close(intervals, max_gap=MAX_ARRAY_GAP):
    if not intervals:
        return []
    out = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s - out[-1][1] <= max_gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def telo_coverage_with_mm(seq, variants, max_mm=1):
    """Fraction of bases inside any variant occurrence allowing <= max_mm mismatches."""
    L = len(seq)
    if L == 0:
        return 0.0
    covered = bytearray(L)
    for v in variants:
        k = len(v)
        if k > L:
            continue
        for i in range(L - k + 1):
            mm = 0
            for j in range(k):
                if seq[i + j] != v[j]:
                    mm += 1
                    if mm > max_mm:
                        break
            if mm <= max_mm:
                for j in range(k):
                    covered[i + j] = 1
    return sum(covered) / L


def classify(seq, motif):
    """Returns (architecture, linker_seq, linker_len, linker_start, linker_end)."""
    F = rotations(motif)
    R = rotations(rc(motif))
    fi = merge_close(find_intervals(seq, F))
    ri = merge_close(find_intervals(seq, R))
    fi = [(s, e) for s, e in fi if e - s >= MIN_ARRAY]
    ri = [(s, e) for s, e in ri if e - s >= MIN_ARRAY]

    if not fi and not ri:
        return "Other", "", 0, -1, -1
    if fi and not ri:
        return "GT-F-AG", "", 0, -1, -1
    if ri and not fi:
        return "GT-R-AG", "", 0, -1, -1

    arrs = [(s, e, "F") for s, e in fi] + [(s, e, "R") for s, e in ri]
    arrs.sort()
    first_kind = arrs[0][2]
    last_kind = arrs[-1][2]

    best = None
    for i in range(len(arrs) - 1):
        if arrs[i][2] == arrs[i + 1][2]:
            continue
        gs, ge = arrs[i][1], arrs[i + 1][0]
        gl = ge - gs
        if gl > 0 and (best is None or gl > best[2]):
            best = (gs, ge, gl, arrs[i][2], arrs[i + 1][2])

    if best is None:
        if first_kind == "F" and last_kind == "R":
            return "GT-F-R-AG", "", 0, -1, -1
        return "Other", "", 0, -1, -1

    gs, ge, gl, prev_kind, next_kind = best
    gap_seq = seq[gs:ge]
    telo_frac = telo_coverage_with_mm(gap_seq, F + R, max_mm=1)
    is_linker = gl >= MIN_LINKER and telo_frac < MAX_LINKER_TELO

    if is_linker:
        if prev_kind == "F" and next_kind == "R":
            return "GT-F-linker-R-AG", gap_seq, gl, gs, ge
        return "GT-R-linker-F-AG", gap_seq, gl, gs, ge

    if first_kind == "F" and last_kind == "R":
        return "GT-F-R-AG", "", 0, -1, -1
    return "Other", "", 0, -1, -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--out-loci", required=True)
    ap.add_argument("--out-kmers", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.final, sep="\t")
    # Use a dict keyed by original row index to avoid misalignment when
    # groupby reorders interleaved rows from different genomes.
    results = {}

    def _store(idx, arch, d6, a6, lseq, llen, l6l, l6r):
        results[idx] = dict(architecture=arch, donor_6mer=d6, acceptor_6mer=a6,
                            linker_seq=lseq, linker_len=llen,
                            linker_left_6mer=l6l, linker_right_6mer=l6r)

    for gid, sub in df.groupby("genome_id", sort=False):
        try:
            fa = find_genome_fasta(gid, args.refseq_dir, args.tara_dir)
        except FileNotFoundError:
            for idx in sub.index:
                _store(idx, "Unknown", "", "", "", 0, "", "")
            continue
        seqs = load_fasta(fa)
        for idx, r in sub.iterrows():
            chrom = seqs.get(r.seqid, "")
            genomic = chrom[r.start - 1:r.end]
            spliced = genomic if r.strand == "+" else rc(genomic)
            motif = r.motif if isinstance(r.motif, str) and r.motif else r.terminal_motif
            if not isinstance(motif, str) or not motif or len(spliced) < 6:
                _store(idx, "Unknown", "", "", "", 0, "", "")
                continue
            arch, lseq, llen, lstart, lend = classify(spliced, motif)
            left6 = right6 = ""
            if lstart >= 2 and lend + 2 <= len(spliced) and llen >= 4:
                left6 = spliced[lstart - 2:lstart + 4]
                right6 = spliced[lend - 4:lend + 2]
            _store(idx, arch, spliced[:6], spliced[-6:], lseq, llen, left6, right6)

    for col in ["architecture", "donor_6mer", "acceptor_6mer",
                "linker_seq", "linker_len", "linker_left_6mer", "linker_right_6mer"]:
        df[col] = [results[i][col] for i in df.index]
    df.to_csv(args.out_loci, sep="\t", index=False)

    rows = []
    for (gid, arch), sub in df.groupby(["genome_id", "architecture"]):
        def emit(side_label, col):
            for kmer, c in sub[col].value_counts().items():
                if kmer:
                    rows.append({"genome_id": gid, "architecture": arch,
                                 "side": side_label, "kmer": kmer, "count": int(c)})
        emit("donor", "donor_6mer")
        emit("acceptor", "acceptor_6mer")
        if arch in ("GT-F-linker-R-AG", "GT-R-linker-F-AG"):
            emit("linker_left", "linker_left_6mer")
            emit("linker_right", "linker_right_6mer")
    pd.DataFrame(rows).to_csv(args.out_kmers, sep="\t", index=False)


if __name__ == "__main__":
    main()
