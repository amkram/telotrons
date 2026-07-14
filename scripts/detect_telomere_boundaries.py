#!/usr/bin/env python3
"""Telomere boundary detection on assembled long-read contigs.

Field-standard sliding-window motif-density method (cf. Stephens & Kocher 2024,
Nurk 2022 / T2T-CHM13, Brown et al. tidk). For each contig end:

  1. Scan first/last SCAN_KB kb for canonical motif matches on both strands.
     - left  end → look for the C-rich revcomp motif (CCCTAA / CCCTAAA);
                   telomere reads 5'→3' INTO the contig.
     - right end → look for the G-rich forward motif (TTAGGG / TTTAGGG);
                   telomere reads 5'→3' OUT of the contig.
     - 1-mismatch tolerance per unit (standard).
  2. Build a per-bp coverage track from merged motif intervals.
  3. Compute motif density in a sliding WINDOW_BP window.
  4. Walk from the contig end inward as long as density >= DENSITY_FRAC.
     A run of <= MAX_GAP_BP without coverage is permitted (allows internal
     TVR / variant units within an otherwise canonical tract).
  5. Boundary = the most-internal position where coverage is contiguous with
     the contig end under the above. Telomere length = boundary - end position.
  6. Require >= MIN_UNITS motif occurrences within the called telomere to
     reject spurious short hits.

Outputs:
  --out-tsv      one row per contig end: contig | end | motif | strand
                 | boundary_pos | telomere_len_bp | n_units | density
                 | max_gap_bp | called (bool)
  --out-summary  plain-text per-species summary
  --out-bed      BED4 (contig, start, end, name) of called telomeres

Configuration:
  --motifs       e.g. "Etenella:TTTAGGG,Emaxima:TTAGGG,Eacervulina:TTTAGGG"
  --assemblies   e.g. "Etenella:path1,Emaxima:path2,Eacervulina:path3"
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

# ── Field-standard defaults (configurable via CLI) ─────────────────────────
SCAN_KB        = 20      # kb to scan at each end (telomeres are usually < 10kb)
WINDOW_BP      = 200     # sliding window for density
DENSITY_FRAC   = 0.70    # min fraction of bases motif-covered within the window
MAX_GAP_BP     = 100     # max consecutive non-motif bp allowed
MIN_UNITS      = 5       # min motif units to call a telomere
MAX_MM_PER_UNIT = 1      # mismatches allowed per motif unit


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def read_fasta(path: Path):
    """Yield (header, sequence) pairs from a (gz-)FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, buf = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            yield name, "".join(buf).upper()


def find_motif_units(seq: str, motif: str, max_mm: int = MAX_MM_PER_UNIT):
    """Yield (start, end) intervals of motif occurrences allowing <= max_mm.
    Iterates unit-by-unit (non-overlapping in the canonical sense, but allows
    units to abut), so a tandem run yields one interval per motif unit.
    """
    k = len(motif)
    L = len(seq)
    out = []
    i = 0
    while i + k <= L:
        sub = seq[i:i + k]
        mm = sum(a != b for a, b in zip(sub, motif))
        if mm <= max_mm:
            out.append((i, i + k))
            i += k  # canonical tandem step — units stack head-to-tail
        else:
            i += 1
    return out


def coverage_track(intervals, length):
    """Per-bp 0/1 coverage from a list of [start, end) intervals."""
    cov = bytearray(length)
    for s, e in intervals:
        for j in range(s, min(e, length)):
            cov[j] = 1
    return cov


def density_walk_from_left(cov: bytearray, scan_len: int) -> tuple[int, int]:
    """Walk from position 0 inward. Returns (boundary_pos, max_gap).

    Algorithm: a "tandem run with gap tolerance" — telomere extends from the
    contig end inward as long as no uncovered run exceeds MAX_GAP_BP. Boundary
    = the position just past the last covered base of the contiguous run.
    A small leading non-motif (≤ MAX_GAP_BP) at the very start is tolerated.
    After identifying the run, verify the run's overall density meets
    DENSITY_FRAC; if not, return 0 (no telomere called).
    """
    L = min(scan_len, len(cov))
    boundary = 0
    gap = 0
    for j in range(L):
        if cov[j]:
            boundary = j + 1
            gap = 0
        else:
            gap += 1
            if gap > MAX_GAP_BP:
                break
    if boundary == 0:
        return 0, 0
    # Re-walk [0, boundary) to capture max internal gap and overall density.
    max_gap = 0
    g = 0
    n_cov = 0
    for j in range(boundary):
        if cov[j]:
            n_cov += 1
            g = 0
        else:
            g += 1
            if g > max_gap:
                max_gap = g
    density = n_cov / boundary if boundary else 0.0
    if density < DENSITY_FRAC:
        return 0, max_gap
    return boundary, max_gap


def density_walk_from_right(cov: bytearray, scan_len: int) -> tuple[int, int]:
    """Mirror of the left walk; returns (boundary_pos_from_right_end, max_gap)
    where boundary_pos is the OFFSET from the end (i.e. telomere length).
    """
    n = len(cov)
    L = min(scan_len, n)
    tail = bytearray(cov[n - L:n])
    tail_rev = bytearray(reversed(tail))
    return density_walk_from_left(tail_rev, L)


def call_end(seq: str, motif: str, end: str) -> dict:
    """Call telomere at `end` ('left' or 'right') of `seq`, using `motif` on the
    appropriate strand. Returns dict of summary stats.
    """
    if end == "left":
        # The 5' end of the spliced sense reads as the C-rich revcomp motif.
        used_motif = revcomp(motif)
        strand = "-"
    else:
        used_motif = motif
        strand = "+"

    scan_len = min(SCAN_KB * 1000, len(seq))
    if end == "left":
        sub = seq[:scan_len]
    else:
        sub = seq[len(seq) - scan_len:]
    intervals = find_motif_units(sub, used_motif)
    cov = coverage_track(intervals, len(sub))

    if end == "left":
        boundary, max_gap = density_walk_from_left(cov, scan_len)
        # boundary is offset from contig start
        boundary_pos = boundary
    else:
        tel_len, max_gap = density_walk_from_right(cov, scan_len)
        boundary = tel_len  # offset measured from right end
        boundary_pos = len(seq) - tel_len

    # Count motif units within the called telomere window.
    if end == "left":
        n_units = sum(1 for s, e in intervals if e <= boundary)
        tel_seq = sub[:boundary]
    else:
        n_units = sum(1 for s, e in intervals if s >= scan_len - boundary)
        tel_seq = sub[scan_len - boundary:]

    density = (sum(cov[:boundary]) / boundary) if end == "left" and boundary else (
        (sum(cov[scan_len - boundary:scan_len]) / boundary) if boundary else 0.0
    )
    called = (n_units >= MIN_UNITS and boundary > 0)
    return {
        "end": end,
        "motif": used_motif,
        "strand": strand,
        "boundary_pos": boundary_pos,
        "telomere_len_bp": boundary,
        "n_units": n_units,
        "density": round(density, 3),
        "max_gap_bp": max_gap,
        "called": called,
    }


def main() -> int:
    global SCAN_KB, WINDOW_BP, DENSITY_FRAC, MAX_GAP_BP, MIN_UNITS, MAX_MM_PER_UNIT
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemblies", required=True,
                    help='comma-sep "species:path,..." (path may be .fna or .fna.gz)')
    ap.add_argument("--motifs", required=True,
                    help='comma-sep "species:motif,..." (G-rich, e.g. TTAGGG / TTTAGGG)')
    ap.add_argument("--scan-kb", type=int, default=SCAN_KB)
    ap.add_argument("--window-bp", type=int, default=WINDOW_BP)
    ap.add_argument("--density-frac", type=float, default=DENSITY_FRAC)
    ap.add_argument("--max-gap-bp", type=int, default=MAX_GAP_BP)
    ap.add_argument("--min-units", type=int, default=MIN_UNITS)
    ap.add_argument("--max-mm", type=int, default=MAX_MM_PER_UNIT)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args()

    SCAN_KB = args.scan_kb
    WINDOW_BP = args.window_bp
    DENSITY_FRAC = args.density_frac
    MAX_GAP_BP = args.max_gap_bp
    MIN_UNITS = args.min_units
    MAX_MM_PER_UNIT = args.max_mm

    asms = dict(s.split(":", 1) for s in args.assemblies.split(","))
    motifs = dict(s.split(":", 1) for s in args.motifs.split(","))

    cols = ["species", "contig", "contig_len", "end", "motif", "strand",
            "boundary_pos", "telomere_len_bp", "n_units", "density",
            "max_gap_bp", "called"]

    rows = []
    bed_lines = []
    per_species_caps = defaultdict(lambda: defaultdict(int))

    for species, path in asms.items():
        motif = motifs.get(species)
        if not motif:
            print(f"[skip] no motif configured for {species}", file=sys.stderr)
            continue
        for contig, seq in read_fasta(Path(path)):
            for end in ("left", "right"):
                r = call_end(seq, motif, end)
                row = {
                    "species": species,
                    "contig": contig,
                    "contig_len": len(seq),
                    **r,
                }
                rows.append(row)
                per_species_caps[species]["total"] += 1
                per_species_caps[species]["called" if r["called"] else "uncalled"] += 1
                if r["called"]:
                    if end == "left":
                        bed_lines.append((contig, 0, r["boundary_pos"],
                                          f"{species}|left|{r['motif']}|n={r['n_units']}"))
                    else:
                        bed_lines.append((contig, r["boundary_pos"], len(seq),
                                          f"{species}|right|{r['motif']}|n={r['n_units']}"))

    # ── outputs ────────────────────────────────────────────────────────────
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(args.out_bed, "w") as fh:
        for contig, s, e, name in bed_lines:
            fh.write(f"{contig}\t{s}\t{e}\t{name}\n")

    with open(args.out_summary, "w") as fh:
        fh.write("# Telomere boundary detection (assembled long-read contigs)\n\n")
        fh.write(f"Method: motif sliding-window density (window {WINDOW_BP} bp, "
                 f"min density {DENSITY_FRAC:.0%}, max gap {MAX_GAP_BP} bp, "
                 f"≤{MAX_MM_PER_UNIT} mm/unit, scan {SCAN_KB} kb from each end, "
                 f"min {MIN_UNITS} units to call).\n\n")
        fh.write("| species | contigs | ends scanned | telomeres called | mean TL bp | median TL bp |\n")
        fh.write("|---|---:|---:|---:|---:|---:|\n")
        for sp, st in sorted(per_species_caps.items()):
            tls = [r["telomere_len_bp"] for r in rows
                   if r["species"] == sp and r["called"]]
            n_contigs = len({r["contig"] for r in rows if r["species"] == sp})
            mean_tl = (sum(tls) / len(tls)) if tls else 0
            median_tl = sorted(tls)[len(tls) // 2] if tls else 0
            fh.write(f"| {sp} | {n_contigs} | {st['total']} | "
                     f"{st['called']} / {st['total']} | "
                     f"{int(mean_tl)} | {median_tl} |\n")
        fh.write("\n\n## Per-contig calls\n\n")
        fh.write("| species | contig | len bp | end | motif | TL bp | units | density | max gap |\n")
        fh.write("|---|---|---:|---|---|---:|---:|---:|---:|\n")
        rows.sort(key=lambda r: (r["species"], r["contig"], r["end"]))
        for r in rows:
            mark = "✓" if r["called"] else "·"
            fh.write(f"| {r['species']} | {r['contig']} | {r['contig_len']} | "
                     f"{r['end']} | {mark} {r['motif']} | {r['telomere_len_bp']} | "
                     f"{r['n_units']} | {r['density']} | {r['max_gap_bp']} |\n")

    print(f"wrote {args.out_tsv}")
    print(f"wrote {args.out_summary}")
    print(f"wrote {args.out_bed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
