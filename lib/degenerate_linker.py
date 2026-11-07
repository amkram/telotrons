"""Anchored, phase-preserving exclusion of degenerate repeats from linker calls."""

from _common import rc
from repeat_policy import coverage, motif_units, rotations

EXTENSION_MISMATCHES = 2
MIN_EXTENSION_UNITS = 3
MIN_CORE_BP = 18


def mismatches(seq, template):
    """Substitutions between two equal-length strings."""
    return sum(base != expected for base, expected in zip(seq, template))


def extend(seq, unit, arm_start, boundary, limit):
    """Furthest end reached past boundary, else boundary.

    Needs a two-unit anchor (<=1 substitution) ending within one unit of boundary, then at least
    MIN_EXTENSION_UNITS same-phase units with <= EXTENSION_MISMATCHES substitutions each.
    """
    width = len(unit)
    best = boundary

    # boundary is a window-union edge, not necessarily a unit boundary,
    # so try anchors ending at boundary and up to width - 1 bp before it
    for anchor_end in range(boundary, max(arm_start + 2 * width - 1, boundary - width), -1):
        anchor_start = anchor_end - 2 * width
        if anchor_start < arm_start:
            continue
        anchor = seq[anchor_start:anchor_end]
        for phase in rotations(unit):
            # anchor: two copies of this phase with at most one substitution
            if set(anchor) - set('ACGT') or mismatches(anchor, phase * 2) > 1:
                continue

            # extend one whole unit at a time, no indels, never past limit
            position = anchor_end
            units = 0
            while position + width <= limit:
                chunk = seq[position : position + width]
                if set(chunk) - set('ACGT') or mismatches(chunk, phase) > EXTENSION_MISMATCHES:
                    break
                units += 1
                position += width

            # keep the anchor and phase that reach furthest
            if units >= MIN_EXTENSION_UNITS and position > best:
                best = position
    return best


def oriented(motif, orientation):
    """The motif as read in an array of the given orientation."""
    return motif if orientation == 'F' else rc(motif)


def evaluate(seq, motif, left, right):
    """Longest non-repeat linker core left after both arrays extend over degenerate repeats."""
    # the linker runs from the left array's end to the right array's start
    linker_start = left['end']
    linker_end = right['start']
    left_reach = extend(
        seq, oriented(motif, left['orientation']), left['start'], linker_start, linker_end
    )

    # the right arm extends leftwards: run on the reversed sequence, then map back
    length = len(seq)
    right_reach = length - extend(
        seq[::-1],
        oriented(motif, right['orientation'])[::-1],
        length - right['end'],
        length - linker_end,
        length - linker_start,
    )

    # residual core measured under the unchanged strict mask
    covered = coverage(seq, motif_units(motif))
    longest = run = 0
    for i in range(left_reach, max(left_reach, right_reach)):
        run = run + 1 if not covered[i] and seq[i] in 'ACGT' else 0
        longest = max(longest, run)
    return dict(
        left_extended_to0=left_reach,
        right_extended_to0=right_reach,
        residual_nonrepeat_bp=longest,
        linker_retained=longest >= MIN_CORE_BP,
    )
