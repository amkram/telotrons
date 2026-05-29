#!/usr/bin/env python3
"""Build the three positive-sequence FASTAs for STREME motif discovery:

  1. telotrons.fa        — concat of every per-genome per-architecture telotron FASTA
  2. non_telo_introns.fa — random sample (--n-non-telo) of introns with
                            telomeric_frac < --non-telo-max-frac, sequences
                            fetched from genome FASTAs via samtools faidx
  3. linkers.fa          — copy of blast_linkers/linker_queries/_all_linkers.fa

Inputs:
  --telotron-fasta-dir   work/results/telotron_fasta/
  --introns-tsv          work/results/all_introns_scanned.tsv  (streamed)
  --linkers-fa           work/results/blast_linkers/linker_queries/_all_linkers.fa
  --refseq-dir, --tara-dir
  --outdir               work/results/streme_inputs/
"""
import argparse
import glob
import os
import random
import shutil
import subprocess

from _common import find_genome_fasta, ensure_uncompressed_faidx


def cat_telotron_fastas(telotron_dir, out_fa):
    n = 0
    with open(out_fa, "w") as oh:
        for fa in sorted(glob.glob(f"{telotron_dir}/*/*.fa")):
            with open(fa) as fh:
                for line in fh:
                    oh.write(line)
                    if line.startswith(">"):
                        n += 1
    return n


def sample_non_telo_introns(introns_tsv, max_frac, n_target, seed):
    """Reservoir sample (gid, seqid, start, end) tuples with telomeric_frac < max_frac."""
    rng = random.Random(seed)
    reservoir = []
    seen = 0
    with open(introns_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx_gid = header.index("genome_id")
        idx_seqid = header.index("seqid")
        idx_start = header.index("start")
        idx_end = header.index("end")
        idx_frac = header.index("telomeric_frac")
        idx_ilen = header.index("intron_len")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                frac = float(f[idx_frac])
                ilen = int(f[idx_ilen])
            except ValueError:
                continue
            if frac >= max_frac or ilen < 20:
                continue
            item = (f[idx_gid], f[idx_seqid], int(f[idx_start]), int(f[idx_end]))
            seen += 1
            if len(reservoir) < n_target:
                reservoir.append(item)
            else:
                j = rng.randint(0, seen - 1)
                if j < n_target:
                    reservoir[j] = item
    return reservoir


def write_non_telo_fa(samples, refseq_dir, tara_dir, work_dir, out_fa):
    os.makedirs(work_dir, exist_ok=True)
    # group by genome to minimise faidx open cost
    by_gid = {}
    for s in samples:
        by_gid.setdefault(s[0], []).append(s)
    n_written = 0
    with open(out_fa, "w") as oh:
        for gid, items in by_gid.items():
            src = find_genome_fasta(gid, refseq_dir, tara_dir, required=False)
            if src is None:
                print(f"  no fasta for {gid} — skip {len(items)} introns")
                continue
            fa = ensure_uncompressed_faidx(src, work_dir)
            regions = [f"{seqid}:{start+1}-{end}" for (_, seqid, start, end) in items]
            # samtools faidx accepts many regions on the command line
            for i in range(0, len(regions), 200):
                chunk = regions[i:i + 200]
                proc = subprocess.run(
                    ["samtools", "faidx", fa] + chunk,
                    check=True, capture_output=True, text=True,
                )
                # re-header with gid prefix so STREME sees unique names
                for line in proc.stdout.splitlines():
                    if line.startswith(">"):
                        oh.write(f">{gid}|{line[1:]}\n")
                        n_written += 1
                    else:
                        oh.write(line + "\n")
    return n_written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--telotron-fasta-dir", required=True)
    ap.add_argument("--introns-tsv", required=True)
    ap.add_argument("--linkers-fa", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-non-telo", type=int, default=5000)
    ap.add_argument("--non-telo-max-frac", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    work_dir = os.path.join(args.outdir, "_work")

    telo_out = os.path.join(args.outdir, "telotrons.fa")
    n_telo = cat_telotron_fastas(args.telotron_fasta_dir, telo_out)
    print(f"telotrons.fa: {n_telo} sequences")

    link_out = os.path.join(args.outdir, "linkers.fa")
    shutil.copyfile(args.linkers_fa, link_out)
    n_link = sum(1 for line in open(link_out) if line.startswith(">"))
    print(f"linkers.fa: {n_link} sequences (copied from {args.linkers_fa})")

    print(f"sampling {args.n_non_telo} introns with frac < {args.non_telo_max_frac}")
    samples = sample_non_telo_introns(
        args.introns_tsv, args.non_telo_max_frac, args.n_non_telo, args.seed
    )
    nt_out = os.path.join(args.outdir, "non_telo_introns.fa")
    n_nt = write_non_telo_fa(samples, args.refseq_dir, args.tara_dir, work_dir, nt_out)
    print(f"non_telo_introns.fa: {n_nt} sequences")


if __name__ == "__main__":
    main()
