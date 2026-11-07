"""Exact tandem telomeric fraction: bases in exact motif rotations (either strand), short
interior gaps bridged.
"""

from _common import rc
from repeat_policy import rotations

GAP = 3  # bridged interior gap (matches curation colouring)


def _cover(seq, variants):
    """Per-base 0/1 list of bases inside any occurrence of any variant, overlaps included."""
    covered = bytearray(len(seq))
    for variant in variants:
        width = len(variant)
        if not width:
            continue
        hit = seq.find(variant)
        while hit != -1:
            for i in range(hit, hit + width):
                covered[i] = 1
            hit = seq.find(variant, hit + 1)
    return covered


def _fill_gaps(labels, gap=GAP):
    """Label every unlabelled interior run of at most gap bases."""
    length = len(labels)
    i = 0
    while i < length:
        if labels[i] == 0:
            run_end = i
            while run_end < length and labels[run_end] == 0:
                run_end += 1

            # fill it if interior (not at either end) and at most gap long
            if i > 0 and run_end < length and (run_end - i) <= gap:
                for position in range(i, run_end):
                    labels[position] = 1
            i = run_end
        else:
            i += 1
    return labels


def telomeric_labels(seq, motif):
    """Per-base 0/1 telomeric labels on either strand, short interior gaps bridged."""
    upper = seq.upper()
    if not upper or not motif:
        return [0] * len(upper)
    motif = motif.upper()
    forward = _cover(upper, rotations(motif))
    reverse = _cover(upper, rotations(rc(motif)))
    labels = [1 if (fwd or rev) else 0 for fwd, rev in zip(forward, reverse)]
    return _fill_gaps(labels)


def exact_telomeric_fraction(seq, motif):
    """Fraction of bases labelled telomeric."""
    if not seq or not motif:
        return 0.0
    return sum(telomeric_labels(seq, motif)) / len(seq)
