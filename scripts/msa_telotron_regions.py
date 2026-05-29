#!/usr/bin/env python3
"""Per-species, per-architecture MAFFT MSAs of telotron regions.

For every (genome_id, architecture) group emit per-region unaligned FASTAs
and aligned MSAs.

Regions per architecture
    All architectures:
        upstream_50       — last 50 bp of the upstream exon
        downstream_50     — first 50 bp of the downstream exon
    Single-arm (GT-F-AG, GT-R-AG):
        intron            — full spliced intron
    Convergent (GT-F-R-AG):
        arm1              — F-arm (5' segment up to the F→R gap)
        arm2              — R-arm (3' segment from the F→R gap)
    Linker (GT-F-linker-R-AG, GT-R-linker-F-AG):
        arm1              — 5' arm before the linker
        linker            — linker sequence itself
        arm2              — 3' arm after the linker

Output layout (per genome / arch):
    <outdir>/<gid>_<organism>/<arch>/
        <region>.fa             ← unaligned input
        <region>.msa.fa         ← MAFFT alignment

Notes:
    - mafft --auto picks L-INS-i / FFT-NS-2 by size; appropriate for our 50-300 bp
      regions without further tuning.
    - region FASTAs with <2 sequences are skipped for alignment (kept as .fa only).
    - All sequences are in display (spliced) orientation; minus-strand loci are
      reverse-complemented so column alignment is meaningful across loci.
"""
import argparse
import os
import subprocess as sp
import sys

import pandas as pd

from _common import rc, slug as _slug, load_fasta, find_genome_fasta

# We reuse the classifier's array-finding logic so arm splits agree with
# what classify_telotron_architecture.py decided when assigning the category.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify_telotron_architecture import rotations, find_intervals, merge_close, rc, MIN_ARRAY


def best_FR_gap(spliced, motif):
    """For convergent (GT-F-R-AG) loci with no real linker: return (gs, ge)
    coordinates of the longest F→R or R→F gap so we can split into arm1/arm2.
    Returns None if no F-then-R or R-then-F transition is found.
    """
    F = rotations(motif)
    R = rotations(rc(motif))
    fi = [(s, e) for s, e in merge_close(find_intervals(spliced, F)) if e - s >= MIN_ARRAY]
    ri = [(s, e) for s, e in merge_close(find_intervals(spliced, R)) if e - s >= MIN_ARRAY]
    if not (fi and ri):
        return None
    arrs = [(s, e, "F") for s, e in fi] + [(s, e, "R") for s, e in ri]
    arrs.sort()
    best = None
    for i in range(len(arrs) - 1):
        if arrs[i][2] == arrs[i + 1][2]:
            continue
        gs, ge = arrs[i][1], arrs[i + 1][0]
        if best is None or ge - gs > best[1] - best[0]:
            best = (gs, ge)
    return best


def regions_for_locus(row, chrom, flank=50):
    """Return dict {region_name: sequence_str} for one locus row.
    Skips regions whose extraction fails (returns empty dict for region key).
    All seqs are in spliced/display orientation.
    """
    L = len(chrom)
    ls = max(0, row.start - 1 - flank)
    le = row.start - 1
    rs = row.end
    re_ = min(L, row.end + flank)
    genomic = chrom[row.start - 1:row.end]
    if row.strand == "-":
        spliced = rc(genomic)
        up = rc(chrom[rs:re_])
        dn = rc(chrom[ls:le])
    else:
        spliced = genomic
        up = chrom[ls:le]
        dn = chrom[rs:re_]

    out = {"upstream_50": up, "downstream_50": dn}
    arch = row.architecture if isinstance(row.architecture, str) else ""

    if arch in ("GT-F-AG", "GT-R-AG"):
        out["intron"] = spliced
    elif arch == "GT-F-R-AG":
        gap = best_FR_gap(spliced, row.motif)
        if gap is not None:
            gs, ge = gap
            out["arm1"] = spliced[:gs]
            out["arm2"] = spliced[ge:]
    elif arch in ("GT-F-linker-R-AG", "GT-R-linker-F-AG"):
        lseq = row.linker_seq if isinstance(row.linker_seq, str) else ""
        if lseq:
            pos = spliced.find(lseq)
            if pos != -1:
                out["arm1"] = spliced[:pos]
                out["linker"] = lseq
                out["arm2"] = spliced[pos + len(lseq):]
    return out


def write_fasta(records, path):
    with open(path, "w") as fh:
        for header, seq in records:
            if not seq:
                continue
            fh.write(">" + header + "\n" + seq + "\n")


def hpc(seq):
    """Homopolymer-compress: collapse runs of identical bases to one base.
    AAAAGGGCCT -> AGCT. Case-preserving on the first base of each run.
    """
    if not seq:
        return seq
    out = [seq[0]]
    for c in seq[1:]:
        if c.upper() != out[-1].upper():
            out.append(c)
    return "".join(out)


def run_mafft(in_fa, out_msa, threads=1):
    """mafft --auto; skip if input has <2 records. Tolerates mafft errors:
    on failure, write an unaligned copy and return False so downstream steps
    still produce something to display.
    """
    n = sum(1 for _ in open(in_fa) if _.startswith(">"))
    if n < 2:
        return False
    cmd = ["mafft", "--auto", "--quiet", "--thread", str(threads), in_fa]
    try:
        with open(out_msa, "w") as fh:
            sp.run(cmd, check=True, stdout=fh, stderr=sp.PIPE)
        return True
    except sp.CalledProcessError as e:
        msg = e.stderr.decode("utf-8", "replace") if e.stderr else "(no stderr)"
        print(f"  mafft failed on {in_fa}: {msg.strip()[:200]}", file=sys.stderr)
        # leave the failed .msa.fa empty so write_region_aln falls back to raw .fa
        try:
            os.remove(out_msa)
        except FileNotFoundError:
            pass
        return False


def read_msa(path):
    """Return list of (header, seq) preserving file order. Joins wrapped lines."""
    recs = []
    cur_h, buf = None, []
    if not os.path.exists(path):
        return recs
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_h is not None:
                    recs.append((cur_h, "".join(buf)))
                cur_h = line[1:]
                buf = []
            else:
                buf.append(line)
        if cur_h is not None:
            recs.append((cur_h, "".join(buf)))
    return recs


# Region order for the combined view, per architecture
def write_region_aln(arch_dir, region, suffix=""):
    """Per-region clean text alignment: fixed-width L#### IDs in left column,
    one line per sequence (no wrapping, no interleaved id rows). Writes
    <region><suffix>.aln.fa  (custom non-FASTA format; using .fa for IDE convenience).
    """
    msa_path = f"{arch_dir}/{region}{suffix}.msa.fa"
    raw_path = f"{arch_dir}/{region}{suffix}.fa"
    recs = read_msa(msa_path) if os.path.exists(msa_path) else read_msa(raw_path)
    if not recs:
        return None
    width = max(len(s) for _, s in recs)
    id_for = {h: f"L{i+1:04d}" for i, (h, _) in enumerate(recs)}
    id_w = max(len(v) for v in id_for.values())

    out_path = f"{arch_dir}/{region}{suffix}.aln.fa"
    label = "homopolymer-compressed" if suffix == ".hpc" else "raw"
    with open(out_path, "w") as fh:
        fh.write(f"# region: {region}   ({label} alignment)\n")
        fh.write(f"# n: {len(recs)}   alignment columns: {width}\n#\n")
        fh.write("# locus index (id -> full header):\n")
        for h in recs:
            fh.write(f"#   {id_for[h[0]]}  {h[0]}\n")
        fh.write("#\n")
        for h, s in recs:
            s_padded = s.ljust(width, "-")
            fh.write(f"{id_for[h]:<{id_w}}  {s_padded}\n")
    return out_path


def write_combined_view_inner(arch_dir, arch, suffix=""):
    """Internal: build combined.aln.fa or combined.hpc.aln.fa from per-region MSAs."""
    regions = REGION_ORDER.get(arch, ["upstream_50", "downstream_50"])
    region_data, region_widths = {}, {}
    for region in regions:
        msa_path = f"{arch_dir}/{region}{suffix}.msa.fa"
        raw_path = f"{arch_dir}/{region}{suffix}.fa"
        recs = read_msa(msa_path) if os.path.exists(msa_path) else read_msa(raw_path)
        if not recs:
            region_data[region], region_widths[region] = {}, 0
            continue
        region_data[region] = {h: s for h, s in recs}
        region_widths[region] = len(recs[0][1])

    seed = next((r for r in regions if region_data[r]), None)
    if seed is None:
        return None
    headers = list(region_data[seed].keys())
    seen = set(headers)
    for r in regions:
        for h in region_data[r]:
            if h not in seen:
                headers.append(h); seen.add(h)

    id_for = {h: f"L{i+1:04d}" for i, h in enumerate(headers)}
    id_w = max((len(v) for v in id_for.values()), default=5)

    fname = "combined.hpc.aln.fa" if suffix == ".hpc" else "combined.aln.fa"
    out_path = f"{arch_dir}/{fname}"
    label = "homopolymer-compressed" if suffix == ".hpc" else "raw"
    with open(out_path, "w") as fh:
        fh.write(f"# architecture: {arch}   ({label} alignment)\n")
        fh.write(f"# region order (space-separated in alignment block):\n")
        fh.write(f"#   {' | '.join(regions)}\n")
        fh.write(f"# region widths (alignment columns): "
                 f"{', '.join(f'{r}={region_widths[r]}' for r in regions)}\n")
        fh.write(f"# n_loci: {len(headers)}\n#\n")
        fh.write(f"# locus index (id -> full header):\n")
        for h in headers:
            fh.write(f"#   {id_for[h]}  {h}\n")
        fh.write("#\n")
        for h in headers:
            parts = []
            for r in regions:
                w = region_widths[r]
                if w == 0:
                    continue
                seq = region_data[r].get(h)
                parts.append(seq if seq is not None else "-" * w)
            fh.write(f"{id_for[h]:<{id_w}}  " + " ".join(parts) + "\n")
    return out_path


REGION_ORDER = {
    "GT-F-AG":          ["upstream_50", "intron", "downstream_50"],
    "GT-R-AG":          ["upstream_50", "intron", "downstream_50"],
    "GT-F-R-AG":        ["upstream_50", "arm1", "arm2", "downstream_50"],
    "GT-F-linker-R-AG": ["upstream_50", "arm1", "linker", "arm2", "downstream_50"],
    "GT-R-linker-F-AG": ["upstream_50", "arm1", "linker", "arm2", "downstream_50"],
    "Other":            ["upstream_50", "downstream_50"],
    "Unknown":          ["upstream_50", "downstream_50"],
}


def write_combined_view(arch_dir, arch):
    """Write both the raw and homopolymer-compressed combined views: one
    fixed-width alignment block per locus, per-region MSAs joined by a single
    space (exon | donor | arm | linker | arm | acceptor | exon, per arch).
    Each view is built independently by write_combined_view_inner."""
    write_combined_view_inner(arch_dir, arch, suffix="")
    write_combined_view_inner(arch_dir, arch, suffix=".hpc")
    return f"{arch_dir}/combined.aln.fa"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch-tsv", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank", type=int, default=50)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.arch_tsv, sep="\t")

    summary = []
    for gid, sub in df.groupby("genome_id", sort=False):
        try:
            fa = find_genome_fasta(gid, args.refseq_dir, args.tara_dir)
        except FileNotFoundError:
            print(f"skip {gid}: FASTA not found")
            continue
        seqs = load_fasta(fa)
        organism = sub.organism.iloc[0]
        slug = f"{_slug(gid)}_{_slug(organism)}"

        # collect per (arch, region) -> list of (header, seq)
        per_arch_region = {}
        for r in sub.itertuples(index=False):
            chrom = seqs.get(r.seqid, "")
            if not chrom:
                continue
            arch = r.architecture if isinstance(r.architecture, str) else "Unknown"
            regs = regions_for_locus(r, chrom, args.flank)
            header = f"{gid}|{r.seqid}|{r.start}-{r.end}|{r.strand}|{arch}"
            for region, seq in regs.items():
                if not seq:
                    continue
                per_arch_region.setdefault((arch, region), []).append((header, seq))

        arches_built = set()
        for (arch, region), recs in per_arch_region.items():
            arch_dir = f"{args.outdir}/{slug}/{arch}"
            os.makedirs(arch_dir, exist_ok=True)

            # raw
            in_fa = f"{arch_dir}/{region}.fa"
            out_msa = f"{arch_dir}/{region}.msa.fa"
            write_fasta(recs, in_fa)
            aligned = run_mafft(in_fa, out_msa, threads=args.threads)

            # homopolymer-compressed
            hpc_recs = [(h, hpc(s)) for h, s in recs]
            in_hpc = f"{arch_dir}/{region}.hpc.fa"
            out_hpc_msa = f"{arch_dir}/{region}.hpc.msa.fa"
            write_fasta(hpc_recs, in_hpc)
            hpc_aligned = run_mafft(in_hpc, out_hpc_msa, threads=args.threads)

            # per-region cleaned alignment text (raw + hpc)
            write_region_aln(arch_dir, region, suffix="")
            write_region_aln(arch_dir, region, suffix=".hpc")

            summary.append({"genome_id": gid, "organism": organism,
                            "architecture": arch, "region": region,
                            "n": len(recs), "aligned": aligned,
                            "hpc_aligned": hpc_aligned})
            print(f"  {gid} {arch} {region}: n={len(recs)} aligned={aligned} hpc={hpc_aligned}")
            arches_built.add(arch)

        # combined per-arch view, after all per-region MSAs are written
        for arch in arches_built:
            arch_dir = f"{args.outdir}/{slug}/{arch}"
            write_combined_view(arch_dir, arch)

    pd.DataFrame(summary).to_csv(f"{args.outdir}/_msa_summary.tsv",
                                  sep="\t", index=False)
    print(f"\nwrote {args.outdir}/_msa_summary.tsv")


if __name__ == "__main__":
    main()
