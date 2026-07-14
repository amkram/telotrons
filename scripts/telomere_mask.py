"""Telomere-fragment masking helper.

Canonical masker for telomeric tandem runs. Use this BEFORE any G/C, G4,
STREME, logo, or BLAST analysis to keep telomere content from leaking
into downstream metrics.

Four prior project hypotheses were retracted because telomere-fragment
content contaminated derived signals:
- G4 density vs telotron count collapsed after masking
  (memory: eimeria_mechanism_g4_supported_2026-05-14)
- Exonic GC enrichment after masking
- RLFS / G-richness signal at subtelomere
  (memory: eimeria_subtelomere_rlfs_2026-05-15)
- C-rich array donor consensus appeared as rotation artefact
  (memory: linker_origin_synthesis_2026-06-03)

Always run `mask_telomere_fragments(seq)` before any composition,
motif, or alignment-based downstream analysis on telomere-adjacent
sequence.
"""

DEFAULT_MOTIFS = ("TTAGGG", "TTTAGGG", "TTAGG", "TTGGGG")


def _revcomp(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _all_rotations(motif):
    """Return every rotation (forward + revcomp) of a motif as a set."""
    rots = set()
    for i in range(len(motif)):
        r = motif[i:] + motif[:i]
        rots.add(r)
        rots.add(_revcomp(r))
    return rots


def mask_telomere_fragments(seq, motifs=DEFAULT_MOTIFS, min_units=2, mask_char="N"):
    """Mask every tandem run of >= min_units copies of any rotation or
    reverse-complement of any motif. Returns a string of the same length
    with masked positions replaced by ``mask_char``.

    Greedy non-overlapping pass; if two motifs overlap on the same window,
    the longest run wins.

    >>> mask_telomere_fragments("ATGCTTAGGGTTAGGGTTAGGGATGC", min_units=2)
    'ATGCNNNNNNNNNNNNNNNNNNATGC'
    >>> mask_telomere_fragments("ATGCTTAGGGATGC", min_units=2)  # one unit only
    'ATGCTTAGGGATGC'
    >>> mask_telomere_fragments("CCCTAACCCTAACCCTAACCCTAA", min_units=4)  # revcomp
    'NNNNNNNNNNNNNNNNNNNNNNNN'
    >>> mask_telomere_fragments("ATTTAGGGTTTAGGGTTTAGGGAT")  # 7-mer x3
    'ANNNNNNNNNNNNNNNNNNNNNAT'
    """
    seq_upper = seq.upper()
    mask = bytearray(b"0" * len(seq_upper))

    # Build the rotation set once. Sort by length DESC so longer motifs win
    # when both could match at the same start position.
    rotation_sets = [(m, _all_rotations(m)) for m in motifs]
    rotation_sets.sort(key=lambda t: -len(t[0]))

    for motif, rots in rotation_sets:
        L = len(motif)
        if L <= 0:
            continue
        i = 0
        n = len(seq_upper)
        while i <= n - L:
            # Anchor at i: try every rotation.
            if seq_upper[i:i + L] in rots:
                # Count how many tandem L-mers follow (each must be in rots).
                run_len = 1
                j = i + L
                while j <= n - L and seq_upper[j:j + L] in rots:
                    run_len += 1
                    j += L
                if run_len >= min_units:
                    for k in range(i, j):
                        mask[k] = ord("1")
                i = j
            else:
                i += 1

    out = "".join(mask_char if mask[k] == ord("1") else seq[k] for k in range(len(seq)))
    return out


def _selftest():
    """Run the docstring assertions inline so this file can be sanity-checked
    without importing doctest. Returns the count of passing assertions."""
    cases = [
        ("ATGCTTAGGGTTAGGGTTAGGGATGC", DEFAULT_MOTIFS, 2,
         "ATGCNNNNNNNNNNNNNNNNNNATGC"),
        ("ATGCTTAGGGATGC", DEFAULT_MOTIFS, 2,
         "ATGCTTAGGGATGC"),
        ("CCCTAACCCTAACCCTAACCCTAA", DEFAULT_MOTIFS, 4,
         "NNNNNNNNNNNNNNNNNNNNNNNN"),
        ("ATTTAGGGTTTAGGGTTTAGGGAT", DEFAULT_MOTIFS, 2,
         "ANNNNNNNNNNNNNNNNNNNNNAT"),
    ]
    n_ok = 0
    for i, (seq, motifs, mu, expected) in enumerate(cases):
        got = mask_telomere_fragments(seq, motifs=motifs, min_units=mu)
        if got == expected:
            n_ok += 1
        else:
            print(f"FAIL case {i}: got {got!r}, expected {expected!r}")
    print(f"selftest: {n_ok}/{len(cases)} pass")
    return n_ok == len(cases)


if __name__ == "__main__":
    import sys
    sys.exit(0 if _selftest() else 1)
