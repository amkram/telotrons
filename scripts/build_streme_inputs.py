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
from telomere_mask import mask_telomere_fragments


def cat_telotron_fastas(telotron_dir, out_fa, mask_telomere=True):
    """Concatenate telotron FASTAs into out_fa.

    Telotron sequences are >50% telomeric repeat by construction (admission
    rules in filter_final_set.py: min_repeat_frac=0.85 single-array,
    bidir_min_repeat_frac=0.40 + bidir_min_hits=3). Running STREME on raw
    telotron FASTAs rediscovers (TTAGGG)n rotations as its top motifs and
    chases the same rotational space across all 10 motifs.

    Four prior project hypotheses were retracted because telomere-fragment
    content leaked into derived signals (G4, exonic GC, RLFS, GTAAACCCTAAA
    "linker donor"). Mask telomere rotations+RC+1-mm before STREME — see
    scripts/telomere_mask.py and memory entries
    eimeria_mechanism_g4_supported_2026-05-14,
    eimeria_subtelomere_rlfs_2026-05-15,
    linker_origin_synthesis_2026-06-03.
    """
    def _emit(oh, header, buf):
        if header is None:
            return 0
        seq = "".join(buf)
        if mask_telomere:
            seq = mask_telomere_fragments(seq)
        oh.write(header)
        for i in range(0, len(seq), 60):
            oh.write(seq[i:i + 60] + "\n")
        return 1

    n = 0
    with open(out_fa, "w") as oh:
        for fa in sorted(glob.glob(f"{telotron_dir}/*/*.fa")):
            with open(fa) as fh:
                header = None
                buf = []
                for line in fh:
                    if line.startswith(">"):
                        n += _emit(oh, header, buf)
                        header = line
                        buf = []
                    else:
                        buf.append(line.strip())
                n += _emit(oh, header, buf)
    return n


def sample_non_telo_introns(introns_tsv, max_frac, n_target, seed, n_per_genome=None):
    """Reservoir-sample (gid, seqid, start, end) tuples with telomeric_frac < max_frac.

    If ``n_per_genome`` is set, runs a per-genome reservoir of that size
    (stratified sampling so STREME isn't over-weighted by species with many
    introns), and ``n_target`` is ignored. Otherwise falls back to a single
    global reservoir of ``n_target``.
    """
    rng = random.Random(seed)
    with open(introns_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx_gid = header.index("genome_id")
        idx_seqid = header.index("seqid")
        idx_start = header.index("start")
        idx_end = header.index("end")
        idx_frac = header.index("telomeric_frac")
        idx_ilen = header.index("intron_len")

        if n_per_genome is not None and n_per_genome > 0:
            # Per-genome reservoir + per-genome seen count.
            reservoirs = {}
            seen_per = {}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                try:
                    frac = float(f[idx_frac])
                    ilen = int(f[idx_ilen])
                except ValueError:
                    continue
                if frac >= max_frac or ilen < 20:
                    continue
                gid = f[idx_gid]
                item = (gid, f[idx_seqid], int(f[idx_start]), int(f[idx_end]))
                seen_per[gid] = seen_per.get(gid, 0) + 1
                seen = seen_per[gid]
                res = reservoirs.setdefault(gid, [])
                if len(res) < n_per_genome:
                    res.append(item)
                else:
                    j = rng.randint(0, seen - 1)
                    if j < n_per_genome:
                        res[j] = item
            out = []
            for gid in sorted(reservoirs):
                out.extend(reservoirs[gid])
            return out

        # Global reservoir (legacy behaviour).
        reservoir = []
        seen = 0
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
    ap.add_argument("--n-non-telo", type=int, default=5000,
                    help="Global reservoir size (used only when --n-per-genome is 0).")
    ap.add_argument("--n-per-genome", type=int, default=0,
                    help="Per-genome reservoir size for stratified sampling "
                         "(>0 disables global sampling so STREME isn't "
                         "over-weighted by species with many introns).")
    ap.add_argument("--non-telo-max-frac", type=float, default=0.10)
    ap.add_argument("--mask-telomere", action="store_true", default=True,
                    help="(Default on.) Mask telomere rotations+RC before STREME.")
    ap.add_argument("--no-mask-telomere", dest="mask_telomere",
                    action="store_false")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    work_dir = os.path.join(args.outdir, "_work")

    telo_out = os.path.join(args.outdir, "telotrons.fa")
    n_telo = cat_telotron_fastas(
        args.telotron_fasta_dir, telo_out, mask_telomere=args.mask_telomere,
    )
    print(f"telotrons.fa: {n_telo} sequences "
          f"(mask_telomere={args.mask_telomere})")

    link_out = os.path.join(args.outdir, "linkers.fa")
    shutil.copyfile(args.linkers_fa, link_out)
    n_link = sum(1 for line in open(link_out) if line.startswith(">"))
    print(f"linkers.fa: {n_link} sequences (copied from {args.linkers_fa})")

    if args.n_per_genome > 0:
        print(f"sampling up to {args.n_per_genome} introns/genome "
              f"with frac < {args.non_telo_max_frac}")
    else:
        print(f"sampling {args.n_non_telo} introns with frac < {args.non_telo_max_frac}")
    samples = sample_non_telo_introns(
        args.introns_tsv, args.non_telo_max_frac, args.n_non_telo, args.seed,
        n_per_genome=(args.n_per_genome if args.n_per_genome > 0 else None),
    )
    nt_out = os.path.join(args.outdir, "non_telo_introns.fa")
    n_nt = write_non_telo_fa(samples, args.refseq_dir, args.tara_dir, work_dir, nt_out)
    print(f"non_telo_introns.fa: {n_nt} sequences")


if __name__ == "__main__":
    main()
