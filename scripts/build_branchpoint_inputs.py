#!/usr/bin/env python3
"""Build branchpoint-search inputs from results/all_introns_scanned.tsv.

For every canonical GT-AG intron with telomeric_frac < threshold and length
≥ min_intron_len, extract two 80-bp windows in splicing orientation:

  positive: positions -80…-1 of the intron (immediately upstream of the AG)
  negative: positions +100…+179 (mid-intron, same intron) — pairs perfectly
            for STREME --objfun de so genome-wide AT bias and simple repeats
            cancel out.

Introns are partitioned into two clade pools:

  tight  — yeast + apicomplexa  (expect near-invariant YACTAAC branchpoint)
  loose  — everything else      (expect degenerate YNYTRAY / YTRAC)

Each pool × {pos, neg} is then dust-masked (BLAST+ dustmasker) so STREME
doesn't latch onto (TA)n / (CACA)n background. Outputs:

  branchpoint_tight_pos.fa  branchpoint_tight_neg.fa
  branchpoint_loose_pos.fa  branchpoint_loose_neg.fa
"""
import argparse
import os
import random
import subprocess as sp

from _common import rc, find_genome_fasta, ensure_uncompressed_faidx

# Yeast (Saccharomyces) + apicomplexa get the "tight" pool — both are known to
# use a near-invariant YACTAAC-style branchpoint. Everything else (animals,
# plants, MAGs, …) goes into "loose".
TIGHT_PREFIXES = ("GCF_000146045.2",   # Saccharomyces cerevisiae
                  "GCF_000006565.2",   # Toxoplasma gondii
                  "GCF_000499385.1",   # Eimeria necatrix
                  "GCF_000499425.1",   # Eimeria acervulina
                  "GCF_000499545.2",   # Eimeria tenella
                  "GCF_000499605.1",   # Eimeria maxima
                  "GCF_000499745.2")   # Eimeria mitis


def clade_of(gid):
    return "tight" if gid in TIGHT_PREFIXES else "loose"


def faidx_batch(fa, regions):
    """Return dict region_key -> sequence (plus strand)."""
    out = {}
    for i in range(0, len(regions), 400):
        chunk = regions[i:i + 400]
        proc = sp.run(["samtools", "faidx", fa] + chunk,
                      check=True, capture_output=True, text=True)
        cur, buf = None, []
        for line in proc.stdout.splitlines():
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur is not None:
            out[cur] = "".join(buf).upper()
    return out


def filter_introns(introns_tsv, max_frac, min_len):
    """Yield (gid, seqid, start, end, strand) for canonical GT-AG introns."""
    with open(introns_tsv) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        IDX = {c: i for i, c in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                frac = float(f[IDX["telomeric_frac"]])
                ilen = int(f[IDX["intron_len"]])
            except ValueError:
                continue
            if frac >= max_frac or ilen < min_len:
                continue
            if f[IDX["splice_class"]] != "GT-AG":
                continue
            yield (f[IDX["genome_id"]], f[IDX["seqid"]],
                   int(f[IDX["start"]]), int(f[IDX["end"]]),
                   f[IDX["strand"]])


def dustmask(in_fa, out_fa):
    """Hard-mask low-complexity bases (lowercase soft-mask → 'N')."""
    soft = in_fa + ".soft"
    sp.run(["dustmasker", "-in", in_fa, "-outfmt", "fasta", "-out", soft],
           check=True)
    # Convert soft-mask (lowercase) → hard-mask (N) so STREME ignores those bases.
    with open(soft) as fh, open(out_fa, "w") as oh:
        for line in fh:
            if line.startswith(">"):
                oh.write(line)
            else:
                oh.write("".join("N" if c.islower() else c
                                  for c in line.rstrip("\n")) + "\n")
    os.remove(soft)


def build_one_species(gid, recs, outdir, work_dir,
                       pos_window, neg_offset, neg_window, min_intron_len,
                       refseq_dir, tara_dir):
    """Extract and dustmask branchpoint windows for a single genome."""
    src = find_genome_fasta(gid, refseq_dir, tara_dir, required=False)
    if src is None:
        print(f"  {gid}: FASTA missing — skip {len(recs)}", flush=True)
        return 0
    fa = ensure_uncompressed_faidx(src, work_dir)
    regions = [f"{seqid}:{start+1}-{end}" for _, seqid, start, end, _ in recs]
    seqs = faidx_batch(fa, regions)
    raw_pos = os.path.join(outdir, f"branchpoint_{gid}_pos.raw.fa")
    raw_neg = os.path.join(outdir, f"branchpoint_{gid}_neg.raw.fa")
    with open(raw_pos, "w") as oh_pos, open(raw_neg, "w") as oh_neg:
        for (_, seqid, start, end, strand), region in zip(recs, regions):
            seq_plus = seqs.get(region, "")
            if not seq_plus:
                continue
            spliced = seq_plus if strand == "+" else rc(seq_plus)
            if len(spliced) < min_intron_len:
                continue
            pos = spliced[-pos_window:]
            neg = spliced[neg_offset: neg_offset + neg_window]
            if len(pos) != pos_window or len(neg) != neg_window:
                continue
            name = f"{gid}|{seqid}:{start}-{end}|{strand}"
            oh_pos.write(f">{name}\n{pos}\n")
            oh_neg.write(f">{name}\n{neg}\n")
    dustmask(raw_pos, os.path.join(outdir, f"branchpoint_{gid}_pos.fa"))
    dustmask(raw_neg, os.path.join(outdir, f"branchpoint_{gid}_neg.fa"))
    os.remove(raw_pos); os.remove(raw_neg)
    n = sum(1 for line in open(os.path.join(outdir, f"branchpoint_{gid}_pos.fa"))
            if line.startswith(">"))
    print(f"  {gid}: {n} window pairs", flush=True)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--introns", required=True,
                    help="results/all_introns_scanned.tsv")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genome-id",
                    help="if given, build inputs for only this genome")
    ap.add_argument("--max-frac", type=float, default=0.10)
    ap.add_argument("--min-intron-len", type=int, default=200)
    ap.add_argument("--pos-window", type=int, default=80)
    ap.add_argument("--neg-offset", type=int, default=100)
    ap.add_argument("--neg-window", type=int, default=80)
    ap.add_argument("--n-per-species", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    work_dir = os.path.join(args.outdir, "_work")
    os.makedirs(work_dir, exist_ok=True)

    rng = random.Random(args.seed)

    # Collect qualifying introns, grouped by genome_id.
    by_gid = {}
    seen_gid = {}
    for rec in filter_introns(args.introns, args.max_frac, args.min_intron_len):
        gid = rec[0]
        if args.genome_id and gid != args.genome_id:
            continue
        seen_gid[gid] = seen_gid.get(gid, 0) + 1
        pool = by_gid.setdefault(gid, [])
        n = args.n_per_species
        if len(pool) < n:
            pool.append(rec)
        else:
            j = rng.randint(0, seen_gid[gid] - 1)
            if j < n:
                pool[j] = rec

    for gid in sorted(by_gid):
        print(f"{gid}: sampled {len(by_gid[gid])} / {seen_gid[gid]} eligible", flush=True)
        build_one_species(gid, by_gid[gid], args.outdir, work_dir,
                          args.pos_window, args.neg_offset, args.neg_window,
                          args.min_intron_len, args.refseq_dir, args.tara_dir)


if __name__ == "__main__":
    main()
