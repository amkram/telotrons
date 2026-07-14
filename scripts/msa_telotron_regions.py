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
    - Two views per region are written: raw (.msa.fa) and motif-collapsed
      (.hpc.msa.fa). The motif-collapsed view replaces every tandem run of
      >=3 forward telomere units with a single ``A+`` token and every reverse
      tandem run with ``A-`` (review fix #39, P0; Brinda 2022 *iScience*
      PMC9633736; Frith 2011 PMC3327517; Bzikadze & Pevzner 2022).
    - WARNING (review findings 7-9, P0): cross-locus MAFFT on tandem-array
      introns does NOT produce orthologous columns; the upstream_50/
      downstream_50 panels are also not column-comparable across host genes.
      Every .aln.fa header carries a column_orthology_warning; never derive
      per-column statistics from these files. Use ortholog MSAs from
      `telotron_ortholog_align` for per-column claims.
"""
import argparse
import os
import subprocess as sp
import sys

import pandas as pd

from _common import rc, slug as _slug, load_fasta, find_genome_fasta

# We reuse the classifier's array-finding logic so arm splits agree with
# what classify_telotron_architecture.py decided when assigning the category.
# Review fix #88 (2026-06-04): drop the local best_FR_gap re-implementation and
# delegate to classifier.classify() for arm boundaries on GT-F-R-AG loci.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify_telotron_architecture import (
    rotations, find_intervals, merge_close, rc, MIN_ARRAY,
    classify as _classify_arch,
)


def best_FR_gap(spliced, motif):
    """For convergent (GT-F-R-AG) loci with no real linker: return (gs, ge)
    coordinates of the F→R or R→F gap used by `classify_telotron_architecture.classify`
    so arm splits agree with the call that assigned the architecture.

    Review fix #88: re-uses the classifier's gap-selection logic via the public
    classify() entry point. Returns None when no convergent transition is found
    or when the classifier did not emit a usable linker span.
    """
    # classify returns (architecture, linker_seq, linker_len, linker_start, linker_end).
    # For convergent (no real linker) inputs the classifier still emits the
    # longest F↔R transition span as linker_start/linker_end on its convergent
    # fallback path; we accept any returned span with lstart>=0.
    _arch, _lseq, _llen, lstart, lend = _classify_arch(spliced, motif)
    if lstart < 0 or lend <= lstart:
        return None
    return (lstart, lend)


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
    """DEPRECATED. Homopolymer compression on telomeric repeats produces
    nonsense (TTAGGG → TAG; review finding 8, Brinda 2022 *iScience*
    PMC9633736). Use `motif_collapse(seq, motif)` instead. Retained only as
    a no-op shim for any external caller; emits the input unchanged.
    """
    return seq


def motif_collapse(seq, motif, min_units=3):
    """Motif-aware tandem-array collapse for cross-locus MSA inputs.

    Replaces every tandem run of >= `min_units` copies of any rotation of
    `motif` with a single forward token `A+`, and every tandem run of any
    rotation of revcomp(motif) with a single reverse token `A-`. Everything
    else (linker, exon flank, intra-array indel) is left untouched.

    This is the recommended substitute for homopolymer compression on tandem
    telomeric arrays (review fix #39, P0). The collapsed sequence keeps
    column-level orthology between arms / linker / flanks across loci that
    differ only in array copy count, which MAFFT --auto otherwise cannot
    align meaningfully (Frith 2011 PMC3327517; Bzikadze & Pevzner 2022
    bioRxiv 2022.09.15.507041).

    Returns the collapsed string. Empty/None input yields empty string.
    """
    if not seq or not motif:
        return seq or ""
    s = seq.upper()
    F = set(rotations(motif))
    R = set(rotations(rc(motif)))
    L = len(motif)
    n = len(s)
    out = []
    i = 0
    while i <= n - L:
        win = s[i:i + L]
        if win in F or win in R:
            in_F = win in F
            run_units = 1
            j = i + L
            while j <= n - L:
                w2 = s[j:j + L]
                # An array switches polarity only via a linker; require all
                # tandem units to match the same orientation set.
                if in_F and w2 in F:
                    run_units += 1
                    j += L
                elif (not in_F) and w2 in R:
                    run_units += 1
                    j += L
                else:
                    break
            if run_units >= min_units:
                out.append("A+" if in_F else "A-")
                i = j
                continue
        out.append(s[i])
        i += 1
    # Trailing tail (last <L bases that could not seed a window) is preserved.
    if i < n:
        out.append(s[i:])
    return "".join(out)


# Sidecar set tracking (arch_dir, region+suffix) pairs whose MAFFT call errored.
# Consumed by write_region_aln / write_combined_view_inner to stamp the
# fallback raw-padded output (review fix #87 / cluster finding 12).
_MAFFT_FAILED = set()


def run_mafft(in_fa, out_msa, threads=1):
    """mafft --auto; skip if input has <2 records. Tolerates mafft errors:
    on failure, write an unaligned copy and return False so downstream steps
    still produce something to display. Failure is recorded in
    `_MAFFT_FAILED` so the .aln.fa fallback is stamped (review fix #87).
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
        # Record the failure so downstream alignment writers can stamp it.
        _MAFFT_FAILED.add(os.path.abspath(in_fa))
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
    used_msa = os.path.exists(msa_path)
    recs = read_msa(msa_path) if used_msa else read_msa(raw_path)
    if not recs:
        return None
    width = max(len(s) for _, s in recs)
    id_for = {h: f"L{i+1:04d}" for i, (h, _) in enumerate(recs)}
    id_w = max(len(v) for v in id_for.values())

    # Review fix #87: detect MAFFT failure so the caller can tell raw-padded
    # output apart from a true alignment.
    raw_abs = os.path.abspath(raw_path)
    alignment_failed = (not used_msa) and (raw_abs in _MAFFT_FAILED)
    # If MAFFT was not run because n<2 or the .msa.fa was missing for any
    # other reason, also flag (conservative): any output sourced from .fa
    # rather than .msa.fa is right-pad-ragged raw sequence, NOT an MSA.
    used_raw_fallback = not used_msa

    out_path = f"{arch_dir}/{region}{suffix}.aln.fa"
    if suffix == ".hpc":
        label = "motif-aware-collapsed"  # review fix #39: was "homopolymer-compressed"
    else:
        label = "raw"
    with open(out_path, "w") as fh:
        fh.write(f"# region: {region}   ({label} alignment)\n")
        fh.write(f"# n: {len(recs)}   alignment columns: {width}\n")
        if alignment_failed:
            fh.write("# alignment_failed: True; showing right-padded raw\n")
        elif used_raw_fallback:
            fh.write("# alignment_skipped: True (n<2 or msa missing); showing right-padded raw\n")
        # Review finding 7/8 (P0): cross-locus MAFFT on un-collapsed tandem
        # arrays does not produce orthologous columns. Even the .hpc / motif-
        # collapsed view is comparable per-arm only after the tandem-run
        # tokenisation; never report per-column statistics from these files.
        fh.write("# column_orthology_warning: cross-locus MSA of tandem-array "
                 "introns; columns are NOT base-orthologous. Do not derive "
                 "per-column statistics. See review_sequence_extraction_and_msa.md.\n#\n")
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
    region_failed = {}
    for region in regions:
        msa_path = f"{arch_dir}/{region}{suffix}.msa.fa"
        raw_path = f"{arch_dir}/{region}{suffix}.fa"
        used_msa = os.path.exists(msa_path)
        recs = read_msa(msa_path) if used_msa else read_msa(raw_path)
        if not recs:
            region_data[region], region_widths[region] = {}, 0
            region_failed[region] = False
            continue
        # Cluster finding 11 (MEDIUM): assert per-record MSA columns equal.
        # If the .msa.fa is malformed (ragged), pad-right to the widest record
        # and record the inconsistency so the combined view stamps a warning.
        w = max(len(s) for _, s in recs)
        ragged = any(len(s) != w for _, s in recs)
        if ragged:
            recs = [(h, s.ljust(w, "-")) for h, s in recs]
        region_data[region] = {h: s for h, s in recs}
        region_widths[region] = w
        raw_abs = os.path.abspath(raw_path)
        region_failed[region] = ragged or ((not used_msa) and (raw_abs in _MAFFT_FAILED))

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
    label = "motif-aware-collapsed" if suffix == ".hpc" else "raw"
    failed_regions = [r for r, f in region_failed.items() if f]
    with open(out_path, "w") as fh:
        fh.write(f"# architecture: {arch}   ({label} alignment)\n")
        fh.write(f"# region order (space-separated in alignment block):\n")
        fh.write(f"#   {' | '.join(regions)}\n")
        fh.write(f"# region widths (alignment columns): "
                 f"{', '.join(f'{r}={region_widths[r]}' for r in regions)}\n")
        if failed_regions:
            fh.write(f"# alignment_failed_regions: {','.join(failed_regions)}; "
                     f"showing right-padded raw for those regions (review fix #87)\n")
        # Review finding 7/8/9 (P0): cross-locus MAFFT on un-collapsed tandem
        # arrays does not produce orthologous columns; flanks_50 across
        # heterologous host genes are not column-comparable either.
        fh.write("# column_orthology_warning: combined view stitches per-region MAFFT\n"
                 "#   alignments across heterologous loci (and across host genes for the\n"
                 "#   upstream_50/downstream_50 panels). Tandem-array regions only become\n"
                 "#   column-comparable after motif_collapse(); flank panels are NOT\n"
                 "#   column-orthologous across host genes. Do not derive per-column\n"
                 "#   statistics. See review_sequence_extraction_and_msa.md findings 7-9.\n")
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
    global MIN_ARRAY
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch-tsv", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank", type=int, default=50)
    ap.add_argument("--threads", type=int, default=4)
    # Must match classify_telotron_architecture's --min-array (single source of
    # truth: config.yaml architecture.min_array_bp, threshold map #4). The
    # Snakefile passes the same value to both rules.
    ap.add_argument("--min-array", type=int, default=MIN_ARRAY,
                    help="min F or R array length (bp) to count as an arm")
    args = ap.parse_args()
    MIN_ARRAY = args.min_array

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

        # collect per (arch, region) -> list of (header, seq, motif). Motif is
        # carried per-record so the motif-aware tandem-collapse path (review
        # fix #39) can use the locus's own canonical motif (TTAGGG vs TTTAGGG
        # vs TTAGG vs ...).
        per_arch_region = {}
        for r in sub.itertuples(index=False):
            chrom = seqs.get(r.seqid, "")
            if not chrom:
                continue
            arch = r.architecture if isinstance(r.architecture, str) else "Unknown"
            regs = regions_for_locus(r, chrom, args.flank)
            header = f"{gid}|{r.seqid}|{r.start}-{r.end}|{r.strand}|{arch}"
            motif = r.motif if isinstance(r.motif, str) and r.motif else "TTAGGG"
            for region, seq in regs.items():
                if not seq:
                    continue
                per_arch_region.setdefault((arch, region), []).append((header, seq, motif))

        arches_built = set()
        for (arch, region), recs in per_arch_region.items():
            arch_dir = f"{args.outdir}/{slug}/{arch}"
            os.makedirs(arch_dir, exist_ok=True)

            # raw: drop the per-record motif before writing the FASTA
            raw_records = [(h, s) for h, s, _m in recs]
            in_fa = f"{arch_dir}/{region}.fa"
            out_msa = f"{arch_dir}/{region}.msa.fa"
            write_fasta(raw_records, in_fa)
            aligned = run_mafft(in_fa, out_msa, threads=args.threads)

            # motif-aware tandem-array collapse (review fix #39, P0). Only
            # apply collapse to regions that can actually contain tandem
            # telomeric arrays (intron / arm* / linker). Flank windows stay
            # uncollapsed.
            tandem_regions = {"intron", "arm1", "arm2", "linker"}
            if region in tandem_regions:
                hpc_recs = [(h, motif_collapse(s, m)) for h, s, m in recs]
            else:
                hpc_recs = list(raw_records)
            in_hpc = f"{arch_dir}/{region}.hpc.fa"
            out_hpc_msa = f"{arch_dir}/{region}.hpc.msa.fa"
            write_fasta(hpc_recs, in_hpc)
            hpc_aligned = run_mafft(in_hpc, out_hpc_msa, threads=args.threads)

            # per-region cleaned alignment text (raw + motif-collapsed)
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
