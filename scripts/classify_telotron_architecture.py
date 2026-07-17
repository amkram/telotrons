#!/usr/bin/env python3
"""Classify each telotron locus by F/R-array architecture and emit boundary 6-mers.

Architectures (6 explicit + Other):
    GT-F-AG            : pure forward (G-rich) array
    GT-R-AG            : pure reverse (C-rich) array
    GT-F-R-AG          : convergent — F then R, no real linker between
    GT-F-linker-R-AG   : F then non-telomeric linker then R
    GT-R-linker-F-AG   : R then non-telomeric linker then F
    Multi-junction     : ≥2 distinct inter-array gaps (review fix #51:
                         do not silently classify by the longest gap; the
                         shorter linker(s) were previously lost to "Other")
    Other              : anything else

A "real linker" must be:
    - length >= MIN_LINKER (default 8 bp — lowered from 15 per
      telotron_two_linker_classes_2026-06-03 to admit the empirically
      bimodal inter-event class with median ~11 bp; review fix #60)
    - <30% bases covered by F/R motif variants allowing 1 mismatch
This kills degenerate-telomeric "linkers" (the trap that bit us previously).

Output:
    --out-loci : per-locus TSV with architecture, donor_6mer, acceptor_6mer,
                 linker_seq, linker_len, linker_start, linker_end
                 (linker_start/linker_end are 0-based indices INTO the
                 spliced-orientation intron; consumed by extract_telotron_fasta
                 and msa_telotron_regions per review_sequence_extraction Fix 4.)
    --out-kmers: long-form TSV (genome_id, architecture, side, kmer, count)
"""
import argparse

import pandas as pd

from _common import rc, rotations, load_fasta, find_genome_fasta

MIN_ARRAY = 18           # minimum length of an F or R array to count
MAX_ARRAY_GAP = 2        # gap (bp) below which two motif hits merge into one array
MIN_LINKER = 8           # minimum gap length to call a linker. Lowered from 15
                         # to admit the inter-event linker class (median 11 bp,
                         # telotron_two_linker_classes_2026-06-03). Central
                         # A−|A+ linkers are ~58 bp and still pass.
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


def telo_coverage_bitvector(seq, variants, max_mm=1):
    """Per-base bitvector (1=inside ≥1 variant with ≤max_mm mm). Review fix #61:
    compute ONCE per locus, then slice for per-window coverage instead of
    re-running the O(L·|variants|·k) scan per linker query."""
    L = len(seq)
    covered = bytearray(L)
    if L == 0:
        return covered
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
    return covered


def telo_coverage_with_mm(seq, variants, max_mm=1):
    """Fraction of bases inside any variant occurrence allowing <= max_mm mismatches.
    Thin wrapper around telo_coverage_bitvector — retained for callers that do
    not pre-compute the locus bitvector."""
    L = len(seq)
    if L == 0:
        return 0.0
    return sum(telo_coverage_bitvector(seq, variants, max_mm)) / L


def coverage_frac_from_bv(bv, start, end):
    """Coverage fraction over bv[start:end] given a pre-computed locus bitvector."""
    if end <= start:
        return 0.0
    return sum(bv[start:end]) / (end - start)


def arm_extents(seq, motif):
    """Return the list of (start, end) F/R arm intervals (>= MIN_ARRAY) in the
    spliced (gene-oriented) intron, merged across F and R and sorted by start.
    Same arm definition `classify()` uses; factored out so telotron_class can be
    computed from arm positions without re-deriving the architecture."""
    F = rotations(motif)
    R = rotations(rc(motif))
    fi = merge_close(find_intervals(seq, F))
    ri = merge_close(find_intervals(seq, R))
    # Motif-scaled minimum: at least MIN_ARRAY bp AND at least 3 motif units.
    # For long-motif taxa (Allium 12bp, Bombus 11bp) the flat 18 bp bar is
    # ~1.5 units and admits scattered noise as arms. 3 units is the accepted
    # minimum for a real telomeric block.
    min_arr = max(MIN_ARRAY, 3 * len(motif))
    arrs = [(s, e) for s, e in fi if e - s >= min_arr] + \
           [(s, e) for s, e in ri if e - s >= min_arr]
    arrs.sort()
    return arrs


def telotron_class(seq, motif, edge=6):
    """Position class of the telomeric array(s) relative to the gene-oriented
    intron boundaries (review-additive column, 2026-06-08):
        complete   : some arm is DONOR-adjacent (5' start <= `edge` bp from the
                     intron start) AND some arm is ACCEPTOR-adjacent (3' end >=
                     intron_len - `edge`). The array spans both splice sites.
        partial    : exactly one side (donor xor acceptor) is reached.
        passenger  : neither boundary is reached — an internal ITS sitting inside
                     a longer, mostly-ancestral intron.
    Returns "none" when no arm >= MIN_ARRAY is found (degenerate / Unknown loci)."""
    L = len(seq)
    arms = arm_extents(seq, motif)
    if not arms:
        return "none"
    donor_adj = any(s <= edge for s, _ in arms)
    acceptor_adj = any(e >= L - edge for _, e in arms)
    if donor_adj and acceptor_adj:
        return "complete"
    if donor_adj or acceptor_adj:
        return "partial"
    return "passenger"


def classify(seq, motif, bv=None):
    """Returns (architecture, linker_seq, linker_len, linker_start, linker_end).

    `bv` is an optional pre-computed F+R 1-mm coverage bitvector for `seq`
    (review fix #61). If omitted, it is computed locally.

    Enumerates ALL inter-array transition gaps (review fix #51): when ≥2
    distinct F↔R gaps pass the (MIN_LINKER, MAX_LINKER_TELO) test, returns
    a `Multi-junction` architecture instead of silently classifying by the
    longest gap and discarding the rest. The reported linker_seq/start/end
    are the LONGEST passing linker (for diagnostics); downstream consumers
    must treat Multi-junction loci as multi-linker.
    """
    F = rotations(motif)
    R = rotations(rc(motif))
    if bv is None:
        bv = telo_coverage_bitvector(seq, F + R, max_mm=1)
    fi = merge_close(find_intervals(seq, F))
    ri = merge_close(find_intervals(seq, R))
    min_arr = max(MIN_ARRAY, 3 * len(motif))
    fi = [(s, e) for s, e in fi if e - s >= min_arr]
    ri = [(s, e) for s, e in ri if e - s >= min_arr]

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

    # Enumerate every F↔R transition gap; classify each as "real linker" or not.
    real_linkers = []   # passes MIN_LINKER + MAX_LINKER_TELO
    all_transitions = []  # any F↔R gap, used for GT-F-R-AG fallback
    for i in range(len(arrs) - 1):
        if arrs[i][2] == arrs[i + 1][2]:
            continue
        gs, ge = arrs[i][1], arrs[i + 1][0]
        gl = ge - gs
        if gl <= 0:
            continue
        prev_kind, next_kind = arrs[i][2], arrs[i + 1][2]
        all_transitions.append((gs, ge, gl, prev_kind, next_kind))
        telo_frac = coverage_frac_from_bv(bv, gs, ge)
        if gl >= MIN_LINKER and telo_frac < MAX_LINKER_TELO:
            real_linkers.append((gs, ge, gl, prev_kind, next_kind))

    # ≥2 real linkers → Multi-junction (review fix #51).
    if len(real_linkers) >= 2:
        gs, ge, gl, _, _ = max(real_linkers, key=lambda x: x[2])
        return "Multi-junction", seq[gs:ge], gl, gs, ge

    if len(real_linkers) == 1:
        gs, ge, gl, prev_kind, next_kind = real_linkers[0]
        gap_seq = seq[gs:ge]
        # A single real linker is only a bidirectional-with-linker call if the
        # locus is exactly two arms (one flanking each side of the linker). If
        # additional arms exist beyond the linker the true structure is a
        # multi-junction, not a linker architecture.
        arm_kinds = [a[2] for a in arrs]
        if len(arm_kinds) == 2:
            if prev_kind == "F" and next_kind == "R":
                return "GT-F-linker-R-AG", gap_seq, gl, gs, ge
            return "GT-R-linker-F-AG", gap_seq, gl, gs, ge
        return "Multi-junction", gap_seq, gl, gs, ge

    # No real linker passed; collapse to a linker-less dyad or Other.
    # BOTH dyad polarities get an explicit label. GT-F-R-AG is the convergent
    # (G-rich 5' -> C-rich 3') dyad; GT-R-F-AG is its divergent mirror image.
    # The mirror case used to fall through to "Other", which hid it from
    # confident_species' BIDIR_ARCHS even though both arms of the telomeric
    # repeat are present — the definition of bidirectional. A genome whose only
    # two telotrons were R-first scored n_bidir=0 and was excluded entirely.
    # The two labels stay DISTINCT rather than merged because dyad orientation
    # is itself a result (~99% of Mixed loci are G-rich-5', memory:
    # mag_tara284_denovo_vs_fusion_2026-06-05); merging them would destroy it.
    if all_transitions and first_kind == "F" and last_kind == "R":
        return "GT-F-R-AG", "", 0, -1, -1
    if all_transitions and first_kind == "R" and last_kind == "F":
        return "GT-R-F-AG", "", 0, -1, -1
    return "Other", "", 0, -1, -1


def main():
    global MAX_LINKER_TELO, MIN_LINKER, MIN_ARRAY
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--out-loci", required=True)
    ap.add_argument("--out-kmers", required=True)
    # Linker definition is threshold-sensitive (and was version-unstable). Expose
    # it so the call is documented/tunable; linker_telomeric_frac is emitted per
    # locus so the F/linker/R classification is auditable, not opaque.
    ap.add_argument("--max-linker-telo", type=float, default=MAX_LINKER_TELO,
                    help="max telomeric fraction (1-mismatch) of a gap to call it a real linker")
    ap.add_argument("--min-linker", type=int, default=MIN_LINKER,
                    help="min gap length (bp) to call a linker")
    # Min array bp must match config.yaml architecture.min_array_bp (single source
    # of truth: top-of-config-file threshold map #4). msa_telotron_regions.py
    # reuses this script's arm definition and accepts the same flag — keep both
    # rules in sync via the Snakefile.
    ap.add_argument("--min-array", type=int, default=MIN_ARRAY,
                    help="min F or R array length (bp) to count as an arm")
    args = ap.parse_args()
    MAX_LINKER_TELO = args.max_linker_telo
    MIN_LINKER = args.min_linker
    MIN_ARRAY = args.min_array

    df = pd.read_csv(args.final, sep="\t")
    # Use a dict keyed by original row index to avoid misalignment when
    # groupby reorders interleaved rows from different genomes.
    results = {}

    # Boundary-6mer naming (review fix #101): the linker is flanked by two
    # array→linker (and linker→array) seams. We label them by directionality
    # along the spliced 5'→3' axis to make the kmer table unambiguous:
    #   array_to_linker_6mer : 2 bp of array + 4 bp of linker (5' seam)
    #   linker_to_array_6mer : 4 bp of linker + 2 bp of array (3' seam)
    # Legacy names (linker_left_6mer / linker_right_6mer) referred to the same
    # bytes but did not state direction; the new names make it explicit.

    def _store(idx, arch, d6, a6, lseq, llen, l6_a2l, l6_l2a, lstart=-1, lend=-1,
               ltf=0.0, tclass="none"):
        results[idx] = dict(architecture=arch, donor_6mer=d6, acceptor_6mer=a6,
                            linker_seq=lseq, linker_len=llen,
                            array_to_linker_6mer=l6_a2l,
                            linker_to_array_6mer=l6_l2a,
                            linker_start=lstart, linker_end=lend,
                            linker_telomeric_frac=ltf,
                            telotron_class=tclass)

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
            # Review fix #61: build the 1-mm F+R coverage bitvector once per
            # locus, then reuse for both classify() and the per-locus
            # linker_telomeric_frac. Saves a full re-scan of the linker.
            variants = rotations(motif) + rotations(rc(motif))
            bv = telo_coverage_bitvector(spliced, variants, max_mm=1)
            arch, lseq, llen, lstart, lend = classify(spliced, motif, bv=bv)
            ltf = coverage_frac_from_bv(bv, lstart, lend) if lseq and lstart >= 0 else 0.0
            tclass = telotron_class(spliced, motif)
            a2l = l2a = ""
            if lstart >= 2 and lend + 2 <= len(spliced) and llen >= 4:
                a2l = spliced[lstart - 2:lstart + 4]
                l2a = spliced[lend - 4:lend + 2]
            _store(idx, arch, spliced[:6], spliced[-6:], lseq, llen, a2l, l2a,
                   lstart, lend, ltf, tclass)

    for col in ["architecture", "donor_6mer", "acceptor_6mer",
                "linker_seq", "linker_len",
                "array_to_linker_6mer", "linker_to_array_6mer",
                "linker_start", "linker_end",
                "linker_telomeric_frac", "telotron_class"]:
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
        if arch in ("GT-F-linker-R-AG", "GT-R-linker-F-AG", "Multi-junction"):
            emit("array_to_linker", "array_to_linker_6mer")
            emit("linker_to_array", "linker_to_array_6mer")
    pd.DataFrame(rows).to_csv(args.out_kmers, sep="\t", index=False)


if __name__ == "__main__":
    main()
