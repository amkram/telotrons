"""Maximal whole-unit tandem tracts, <=1 substitution per adjacent unit pair, 0-based half-open."""

from repeat_policy import motif_units


def tracts(seq, motif, min_units=4):
    """(start, end, mismatches, pattern) of each maximal tract, <=1 substitution per unit pair."""
    upper = seq.upper()
    width = len(motif)
    found = []

    def flush(blocks, pattern):
        # emit a tract of >= min_units (start, mismatches) blocks, a lone unit only if exact
        if len(blocks) >= min_units and (len(blocks) > 1 or blocks[0][1] == 0):
            total = sum(unit_mismatches for _, unit_mismatches in blocks)
            found.append((blocks[0][0], blocks[-1][0] + width, total, pattern))

    # scan each rotation/strand pattern in each reading frame
    for pattern in sorted(motif_units(motif)):
        for frame in range(width):
            blocks = []
            for start in range(frame, len(upper) - width + 1, width):
                unit = upper[start : start + width]
                unit_mismatches = sum(base != expected for base, expected in zip(unit, pattern))

                # an ambiguous unit or one with >1 mismatch ends the tract
                if set(unit) - set('ACGT') or unit_mismatches > 1:
                    flush(blocks, pattern)
                    blocks = []
                    continue

                # two adjacent units may hold at most one substitution together
                if blocks and blocks[-1][1] + unit_mismatches > 1:
                    flush(blocks, pattern)
                    blocks = []
                blocks.append((start, unit_mismatches))
            flush(blocks, pattern)
    return sorted(set(found))
