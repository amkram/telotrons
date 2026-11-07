"""Repeat evidence: exact units, or same-phase two-unit windows with <=1 substitution."""

from functools import lru_cache

from _common import rc


@lru_cache(maxsize=256)
def pattern_sets(units):
    """({width: exact units}, {width: two-unit windows with at most one substitution})."""
    exact, paired = {}, {}
    for unit in units:
        unit = unit.upper()
        if not unit or set(unit) - set('ACGT'):
            raise ValueError('Repeat motifs must contain resolved ACGT bases')
        width = len(unit)
        exact.setdefault(width, set()).add(unit)

        # the two-unit window plus all its single-substitution variants
        window = unit * 2
        near = paired.setdefault(2 * width, set())
        near.add(window)
        for i in range(len(window)):
            for base in 'ACGT':
                if base != window[i]:
                    near.add(window[:i] + base + window[i + 1 :])
    return exact, paired


def matching_spans(seq, units):
    """(start, end, is_exact) of every window matching a unit or a near two-unit window."""
    upper = seq.upper()
    # sorted tuple keeps the lru_cache key hashable and order-free
    exact, paired = pattern_sets(tuple(sorted(set(units))))
    for table, is_exact in ((exact, True), (paired, False)):
        for width, patterns in table.items():
            for start in range(len(upper) - width + 1):
                if upper[start : start + width] in patterns:
                    yield start, start + width, is_exact


def coverage(seq, units):
    """Per-base 0/1 mask of bases inside any matching window."""
    covered = bytearray(len(seq))
    for start, end, _ in matching_spans(seq, units):
        covered[start:end] = b'\x01' * (end - start)
    return covered


def rotations(unit):
    """Every cyclic rotation of a unit."""
    return {unit[i:] + unit[:i] for i in range(len(unit))}


def motif_units(motif):
    """Rotations of the motif and of its reverse complement."""
    motif = motif.upper()
    return rotations(motif) | rotations(rc(motif))


def repeat_labels(seq, motifs):
    """Per-base labels: F/R exact unit, f/r near two-unit window, '.' otherwise."""
    labels = ['.'] * len(seq)

    # F collects motif rotations, R reverse-complement rotations
    units_by_tag = {'F': set(), 'R': set()}
    for motif in motifs:
        motif = motif.upper()
        if not motif:
            continue
        units_by_tag['F'].update(rotations(motif))
        units_by_tag['R'].update(rotations(rc(motif)))
    spans = {tag: list(matching_spans(seq, units)) for tag, units in units_by_tag.items()}

    # precedence exact > approximate, then F > R (figure convention)
    # a base keeps the first label it receives
    for want_exact in (True, False):
        for tag in ('F', 'R'):
            mark = tag if want_exact else tag.lower()
            for start, end, is_exact in spans[tag]:
                if is_exact == want_exact:
                    for i in range(start, end):
                        if labels[i] == '.':
                            labels[i] = mark
    return ''.join(labels)
