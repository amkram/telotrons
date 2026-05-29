#!/usr/bin/env python3
"""Per-genome BED of unannotated ORFs to exclude from the interstitial array set.

Six-frame ATG→stop scan: any ORF >= --min-orf-nt (default 450 nt = 150 codons)
is added to the exclusion BED. On a shuffled E. necatrix genome, ORFs >= 300 nt
cover 10.6% of bp (indistinguishable from noise); >= 450 nt cover 1.9% on the
shuffled control vs 5.3% on the real genome, so 450 is the noise floor.

Per genome we emit `{outdir}/{genome_id}.bed` (sorted, merged). Empty file if
the FASTA is missing.
"""
import argparse
import gzip
import os
import subprocess as sp
import sys
import tempfile

import pandas as pd

from _common import find_genome_fasta, maybe_decompress, ensure_tool_on_path

# Make required binaries reachable for bare `python scripts/...` runs outside
# `snakemake --use-conda` (no-op when they're already on PATH).
ensure_tool_on_path("bedtools", "samtools")


COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
STOPS = frozenset((b"TAA", b"TAG", b"TGA"))


def rc_bytes(b):
    return b.translate(COMP)[::-1]


def iter_fasta_bytes(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        sid, buf = None, []
        for raw in fh:
            line = raw.rstrip()
            if not line:
                continue
            if line[:1] == b">":
                if sid is not None:
                    yield sid, b"".join(buf).upper()
                sid = line[1:].split(None, 1)[0].decode()
                buf = []
            else:
                buf.append(line)
        if sid is not None:
            yield sid, b"".join(buf).upper()


def find_orfs_one_frame(seq, frame, min_nt):
    out = []
    L = len(seq)
    atg = -1
    i = frame
    while i + 3 <= L:
        codon = seq[i:i + 3]
        if codon in STOPS:
            if atg >= 0 and (i + 3 - atg) >= min_nt:
                out.append((atg, i + 3))
            atg = -1
        elif atg < 0 and codon == b"ATG":
            atg = i
        i += 3
    return out


def orf_bed_for_contig(chrom, seq, min_nt):
    L = len(seq)
    rows = []
    for frame in range(3):
        for s, e in find_orfs_one_frame(seq, frame, min_nt):
            rows.append(f"{chrom}\t{s}\t{e}\t.\t.\t+")
    rseq = rc_bytes(seq)
    for frame in range(3):
        for s, e in find_orfs_one_frame(rseq, frame, min_nt):
            rows.append(f"{chrom}\t{L - e}\t{L - s}\t.\t.\t-")
    return rows


def orf_bed(fa_path, min_nt):
    rows = []
    for sid, seq in iter_fasta_bytes(fa_path):
        rows.extend(orf_bed_for_contig(sid, seq, min_nt))
    return "\n".join(rows) + ("\n" if rows else "")


def sort_and_merge(bed_text, fai_path):
    if not bed_text.strip():
        return ""
    sort = sp.run(["bedtools", "sort", "-faidx", fai_path, "-i", "-"],
                  input=bed_text, capture_output=True, text=True, check=True)
    merged = sp.run(["bedtools", "merge", "-i", "-"],
                    input=sort.stdout, capture_output=True, text=True, check=True)
    return merged.stdout


def process_genome(gid, fa_path, min_orf, outdir):
    out_path = os.path.join(outdir, f"{gid}.bed")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        print(f"  {gid}: mask already present, skipping", flush=True)
        return
    with tempfile.TemporaryDirectory(prefix=f"mask_{gid}_") as tmpdir:
        fa_plain = maybe_decompress(fa_path, tmpdir)
        sp.run(["samtools", "faidx", fa_plain], check=True,
               stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        fai = fa_plain + ".fai"
        of_bed = orf_bed(fa_plain, min_orf)
        merged = sort_and_merge(of_bed, fai)
        with open(out_path, "w") as fh:
            fh.write(merged)
        n_orfs = of_bed.count("\n")
        n_merged = merged.count("\n")
        print(f"  {gid}: {n_orfs} ORFs → {n_merged} merged intervals", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-orf-nt", type=int, default=450)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t")
    for _, m in manifest.iterrows():
        gid = m.genome_id
        fa = find_genome_fasta(gid, args.refseq_dir, args.tara_dir, required=False)
        if not fa:
            print(f"skip {gid}: no FASTA", flush=True)
            open(os.path.join(args.outdir, f"{gid}.bed"), "w").close()
            continue
        try:
            process_genome(gid, fa, args.min_orf_nt, args.outdir)
        except sp.CalledProcessError as e:
            print(f"  {gid}: tool failed ({e.cmd[0]} exit {e.returncode})", flush=True)


if __name__ == "__main__":
    sys.exit(main())
