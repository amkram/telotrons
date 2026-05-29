#!/usr/bin/env python3
"""Find interstitial (non-genic, non-intronic, non-terminal) telomeric arrays.

Pipeline per genome
  1. seqkit locate -p MOTIFS --bed   → all motif hits on both strands
  2. bedtools sort | bedtools merge -d MAX_GAP -s   → strand-aware arrays
  3. filter arrays by MIN_ARRAY_LEN
  4. flag terminal arrays (within TERMINAL_BP of contig end) and drop them
  5. derive intron BED from the GFF via `gt gff3 -addintrons`
  6. extract gene BED from the GFF (awk over $3=="gene")
  7. bedtools intersect -v arrays vs (genes ∪ introns)   → non-genic, non-intronic
  8. emit each survivor's first6/first12/last6/last12 (plus strand 5'→3')

Tools used: seqkit locate (motif scan), bedtools (interval ops),
GenomeTools `gt gff3 -addintrons` (intron derivation) — same toolchain the
main pipeline relies on for telotron scanning.

Output TSV columns:
    genome_id organism group seqid start end strand array_len
    first6 first12 last6 last12 distance_to_contig_end contig_len
    motif cat5_3 upstream_flank10 first_unit last_unit downstream_flank10

`cat5_3` is one of {GG, GC, CG, CC, NG, GN, …}: the orientation of the 5' and
3' ends of the array (G-rich-telomeric vs C-rich-telomeric), classified
independently. Hybrid arrays (GC, CG) come from strand-switched repeats.
"""
import argparse
import gzip
import os
import subprocess as sp
import sys
import tempfile

import pandas as pd

from _common import (rc, expand_motifs, find_genome_files, maybe_decompress,
                     load_fasta, ensure_tool_on_path)

# Make required binaries reachable for bare `python scripts/...` runs outside
# `snakemake --use-conda` (no-op when they're already on PATH).
ensure_tool_on_path("bedtools", "seqkit", "samtools")

MAX_GAP = 12           # merge motif hits within 12 bp (≈2 motif units) into one array
MIN_ARRAY_LEN = 30     # minimum array length to consider
TERMINAL_BP = 5000     # arrays within 5 kb of a contig end are "terminal", excluded


def contig_lengths(fai_path):
    out = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            out[parts[0]] = int(parts[1])
    return out


def seqkit_locate_bed(fa_path, motifs_csv, threads):
    """Return BED text of motif hits on the PLUS strand only.

    Our expanded motif list already includes reverse complements, so scanning
    both strands would just double-count every array.
    """
    cmd = ["seqkit", "locate", "--bed", "--ignore-case", "--only-positive-strand",
           "--threads", str(threads),
           "-p", motifs_csv, fa_path]
    res = sp.run(cmd, check=True, capture_output=True, text=True)
    return res.stdout


def bedtools_merge(bed_text, fai_path, max_gap):
    """Sort + strand-UNaware merge (we've collapsed to plus-strand hits already)."""
    sort = sp.run(["bedtools", "sort", "-faidx", fai_path, "-i", "-"],
                  input=bed_text, check=True, capture_output=True, text=True)
    merged = sp.run(["bedtools", "merge", "-d", str(max_gap), "-i", "-"],
                    input=sort.stdout, check=True, capture_output=True, text=True)
    return merged.stdout


def gene_bed_from_gff(gff_path):
    """BED of every record with feature == gene (any source)."""
    out = []
    opener = gzip.open if gff_path.endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            try:
                s, e = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            out.append(f"{f[0]}\t{s}\t{e}\t.\t.\t{f[6]}")
    return "\n".join(out) + ("\n" if out else "")


def intron_bed_via_gt(gff_path, tmpdir):
    """Use gt gff3 -addintrons to derive intron BED."""
    addintrons = os.path.join(tmpdir, "addintrons.gff")
    with open(addintrons, "w") as fh:
        sp.run(["gt", "gff3", "-addintrons", "-tidy", "-retainids", gff_path],
               check=True, stdout=fh, stderr=sp.DEVNULL)
    out = []
    with open(addintrons) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "intron":
                continue
            try:
                s, e = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            out.append(f"{f[0]}\t{s}\t{e}\t.\t.\t{f[6]}")
    return "\n".join(out) + ("\n" if out else "")


def bedtools_intersect_v(array_bed, excl_bed, tmpdir):
    if not excl_bed.strip():
        return array_bed
    excl_path = os.path.join(tmpdir, "exclude.bed")
    with open(excl_path, "w") as fh:
        fh.write(excl_bed)
    res = sp.run(["bedtools", "intersect", "-v", "-a", "-", "-b", excl_path],
                 input=array_bed, check=True, capture_output=True, text=True)
    return res.stdout


def extend_array(seq, start, end, motif_list):
    """Extend an interstitial array outward, unit-by-unit, allowing ≤1 mismatch per unit.

    For each candidate motif (and its reverse complement), pick the rotation
    that best matches the array's first motif-length window — that fixes the
    tandem-repeat phase. Then walk outward in blocks of len(motif): accept a
    block as part of the array iff it differs from the phased-expected unit by
    at most one base. Return the (start, end) that produced the longest array
    across all motifs tried.
    """
    L_total = len(seq)
    best = (start, end, None)
    best_len = end - start

    candidates = set()
    for m in motif_list:
        mu = m.upper()
        candidates.add(mu)
        candidates.add(rc(mu))

    for m in candidates:
        L = len(m)
        if L == 0 or end - start < L:
            continue
        rots = [m[i:] + m[:i] for i in range(L)]
        # Anchor rotation: best match to seq[start:start+L]
        head = seq[start:start + L]
        anchor = max(rots, key=lambda r: sum(a == b for a, b in zip(head, r)))
        # Bail if even the anchor frame is a poor fit (>1 mismatch on the head).
        head_mm = sum(a != b for a, b in zip(head, anchor))
        if head_mm > 1:
            continue

        # Extend right: expected block at offset k from `start` is
        # anchor[(k + j) % L] for j in 0..L-1; for k = (end - start) this
        # is a rotation of `anchor`.
        e2 = end
        while e2 + L <= L_total:
            k = e2 - start
            expected = anchor[k % L:] + anchor[:k % L]
            block = seq[e2:e2 + L]
            if sum(a != b for a, b in zip(block, expected)) <= 1:
                e2 += L
            else:
                break

        # Extend left: block immediately before `s2` is at offset (s2 - L - start)
        # from the anchor; that's a negative offset whose mod-L is well-defined.
        s2 = start
        while s2 - L >= 0:
            k = (s2 - L - start) % L
            expected = anchor[k:] + anchor[:k]
            block = seq[s2 - L:s2]
            if sum(a != b for a, b in zip(block, expected)) <= 1:
                s2 -= L
            else:
                break

        if e2 - s2 > best_len:
            best = (s2, e2, m)
            best_len = e2 - s2
    return best


def edge_orient(window, motif_list):
    """Classify a sequence window as G-rich-telomeric, C-rich-telomeric, or N.

    Returns 'G' if the window better matches any rotation of a base motif,
    'C' if it better matches its reverse complement, 'N' on a tie / no signal.
    Used to classify the 5' and 3' edges of an interstitial array independently
    so hybrid arrays (strand-switched within the merged interval) can be split
    out as GC or CG, alongside pure GG / CC.
    """
    if not window:
        return "N"
    g_best, c_best = -1, -1
    for m in motif_list:
        L = len(m)
        if L == 0 or len(window) < L:
            continue
        head = window[:L]
        for i in range(L):
            r = m[i:] + m[:i]
            g_best = max(g_best, sum(a == b for a, b in zip(head, r)))
            rcr = rc(r)
            c_best = max(c_best, sum(a == b for a, b in zip(head, rcr)))
    if g_best > c_best:
        return "G"
    if c_best > g_best:
        return "C"
    return "N"


def parse_array_bed(bed_text):
    """Yield dicts: chrom, start, end (strand-unaware after merge)."""
    for line in bed_text.splitlines():
        if not line: continue
        f = line.split("\t")
        if len(f) < 3: continue
        yield {"chrom": f[0], "start": int(f[1]), "end": int(f[2])}


def process_genome(gid, organism, group, fa_path, gff_path, motifs_csv,
                   threads, terminal_bp, min_len, extra_mask_bed=None,
                   base_motifs=None):
    """Returns list[dict] of interstitial array records (or [] if errors)."""
    with tempfile.TemporaryDirectory(prefix="interst_") as tmpdir:
        fa_plain = maybe_decompress(fa_path, tmpdir)
        # samtools faidx for sort/length
        sp.run(["samtools", "faidx", fa_plain], check=True,
               stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        fai = fa_plain + ".fai"
        clen = contig_lengths(fai)

        bed_text = seqkit_locate_bed(fa_plain, motifs_csv, threads)
        if not bed_text.strip():
            return []
        merged = bedtools_merge(bed_text, fai, MAX_GAP)

        # filter by length + terminal
        kept_lines = []
        for r in parse_array_bed(merged):
            L = r["end"] - r["start"]
            if L < min_len:
                continue
            seq_len = clen.get(r["chrom"], 0)
            dist_end = min(r["start"], seq_len - r["end"])
            if dist_end < terminal_bp:
                continue  # terminal — skip
            kept_lines.append(f"{r['chrom']}\t{r['start']}\t{r['end']}\t.\t{L}\t+")
        if not kept_lines:
            return []
        array_bed = "\n".join(kept_lines) + "\n"

        # genic/intronic exclude (annotated) + unannotated-mask (miniprot ∪ ORFs)
        excl = ""
        if gff_path:
            excl += gene_bed_from_gff(gff_path)
            try:
                excl += intron_bed_via_gt(gff_path, tmpdir)
            except (sp.CalledProcessError, FileNotFoundError):
                pass  # GFF may not be gt-clean (or gt missing); gene-only filter still applied
        if extra_mask_bed and os.path.exists(extra_mask_bed) and os.path.getsize(extra_mask_bed) > 0:
            # The mask is 3-col (chrom start end); pad to 6-col so bedtools
            # accepts a uniform exclusion file alongside the gene/intron rows.
            with open(extra_mask_bed) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 6:
                        parts = parts + ["."] * (5 - len(parts)) + ["+"]
                    excl += "\t".join(parts[:6]) + "\n"
        survivors = bedtools_intersect_v(array_bed, excl, tmpdir)

        if not survivors.strip():
            return []

        # extract first/last k-mers; infer orientation (G-rich = +, C-rich = -)
        fa_seqs = load_fasta(fa_plain)
        out = []
        for line in survivors.splitlines():
            if not line: continue
            f = line.split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            seq = fa_seqs.get(chrom, "")
            L = len(seq)
            if L == 0: continue
            # Extend boundaries unit-by-unit with ≤1 mm/unit, picking the
            # rotation/motif that yields the longest array.
            best_motif = None
            if base_motifs:
                start, end, best_motif = extend_array(seq, start, end, base_motifs)
            arr = seq[start:end]
            arr_len = end - start
            dist_end = min(start, L - end)

            # Orientation: if more G+T than A+C, plus strand is the G-rich strand.
            # Otherwise reorient (rc) so first/last k-mers reflect the canonical
            # telomeric (G-rich) orientation.
            gt = arr.count("G") + arr.count("T")
            ac = arr.count("A") + arr.count("C")
            if gt >= ac:
                strand = "+"
                aseq = arr
            else:
                strand = "-"
                aseq = rc(arr)

            # Boundary k-mers: extend 2 bp into the flanking sequence on each side
            # so first6/first12 start 2 bp before the array; last6/last12 end 2 bp after.
            ext_s = max(0, start - 2)
            ext_e = min(L, end + 2)
            arr_ext = seq[ext_s:ext_e]
            if gt >= ac:
                aseq_ext = arr_ext
            else:
                aseq_ext = rc(arr_ext)

            # Sequence-logo windows: 10 bp upstream flank ++ first unit, and
            # last unit ++ 10 bp downstream flank, in the canonical G-rich
            # orientation. Padded with 'N' when the contig edge is closer than
            # 10 bp (rare, since terminal arrays are already excluded).
            motif_len = len(best_motif) if best_motif else 0
            cat5_3 = ""
            if motif_len > 0 and arr_len >= motif_len and base_motifs:
                # Classify each end of the *output* (already-reoriented) array
                # against the configured base motifs. Each end is independent
                # so strand-switched hybrid arrays get GC or CG.
                o5 = edge_orient(aseq[:motif_len], base_motifs)
                o3 = edge_orient(aseq[-motif_len:], base_motifs)
                cat5_3 = f"{o5}{o3}"
            if motif_len > 0 and arr_len >= motif_len:
                up_raw = seq[max(0, start - 10):start]
                dn_raw = seq[end:end + 10]
                up_raw = ("N" * (10 - len(up_raw))) + up_raw
                dn_raw = dn_raw + ("N" * (10 - len(dn_raw)))
                first_unit_pl = seq[start:start + motif_len]
                last_unit_pl  = seq[end - motif_len:end]
                if strand == "+":
                    upstream_flank10   = up_raw
                    downstream_flank10 = dn_raw
                    first_unit = first_unit_pl
                    last_unit  = last_unit_pl
                else:
                    upstream_flank10   = rc(dn_raw)
                    downstream_flank10 = rc(up_raw)
                    first_unit = rc(last_unit_pl)
                    last_unit  = rc(first_unit_pl)
            else:
                upstream_flank10 = downstream_flank10 = first_unit = last_unit = ""

            out.append({
                "genome_id": gid, "organism": organism, "group": group,
                "seqid": chrom, "start": start, "end": end, "strand": strand,
                "array_len": arr_len,
                "first6": aseq_ext[:6], "first12": aseq_ext[:12],
                "last6": aseq_ext[-6:], "last12": aseq_ext[-12:],
                "distance_to_contig_end": dist_end, "contig_len": L,
                "motif": (best_motif or ""),
                "cat5_3": cat5_3,
                "upstream_flank10": upstream_flank10,
                "first_unit": first_unit,
                "last_unit": last_unit,
                "downstream_flank10": downstream_flank10,
            })
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="all_genomes.tsv")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--motifs", required=True, help="comma-separated motifs")
    ap.add_argument("--out", required=True, help="interstitial arrays TSV")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--terminal-bp", type=int, default=TERMINAL_BP)
    ap.add_argument("--min-array-len", type=int, default=MIN_ARRAY_LEN)
    ap.add_argument("--mask-dir", default=None,
                    help="optional dir of per-genome BED masks (miniprot ∪ ORFs); "
                         "files named {genome_id}.bed are added to the exclusion set")
    args = ap.parse_args()

    motif_list = [m.strip() for m in args.motifs.split(",") if m.strip()]
    expanded = expand_motifs(motif_list)
    motifs_csv = ",".join(expanded)

    manifest = pd.read_csv(args.manifest, sep="\t")
    all_rows = []
    for _, m in manifest.iterrows():
        gid, organism, group = m.genome_id, m.organism, m.group
        fa, gff = find_genome_files(gid, args.refseq_dir, args.tara_dir)
        if not fa:
            print(f"skip {gid}: no FASTA", flush=True)
            continue
        extra_mask = None
        if args.mask_dir:
            cand = os.path.join(args.mask_dir, f"{gid}.bed")
            if os.path.exists(cand):
                extra_mask = cand
        try:
            rows = process_genome(gid, organism, group, fa, gff, motifs_csv,
                                  args.threads, args.terminal_bp, args.min_array_len,
                                  extra_mask_bed=extra_mask, base_motifs=motif_list)
        except sp.CalledProcessError as e:
            print(f"  {gid}: tool failed ({e.cmd[0]} exit {e.returncode})", flush=True)
            continue
        tags = []
        if not gff: tags.append("no GFF")
        if not extra_mask: tags.append("no extra mask")
        tag_str = f"({'; '.join(tags)})" if tags else ""
        print(f"  {gid}: {len(rows)} interstitial arrays {tag_str}", flush=True)
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {args.out} ({len(df)} arrays across {df.genome_id.nunique() if len(df) else 0} genomes)")


if __name__ == "__main__":
    sys.exit(main())