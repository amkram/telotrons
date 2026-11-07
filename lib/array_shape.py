"""Telomeric repeat arrays and intron architecture, 0-based half-open.

An array joins same-orientation repeat windows across gaps <= MAX_GAP_BP, gap bases unsupported.
"""

from _common import rc
from repeat_policy import matching_spans, rotations

MAX_GAP_BP = 3
EDGE_BP = 8
MIN_UNITS = 3
MIN_LINKER_BP = 18


def tier(covered_bp, length):
    """'standard' (primary) above 85% coverage, 'additional' above 50% to 85%, else 'excluded'."""
    if length <= 0 or not 0 <= covered_bp <= length:
        raise ValueError('Invalid covered-base count or length')
    if covered_bp * 100 > 85 * length:
        return 'standard'
    if covered_bp * 2 > length:
        return 'additional'
    return 'excluded'


def runs(mask):
    """(start, end) of each nonzero run of a mask."""
    start = None
    # the appended 0 closes a final run
    for i, value in enumerate(list(mask) + [0]):
        if value and start is None:
            start = i
        elif not value and start is not None:
            yield start, i
            start = None


def components(spans, orientation):
    """Spans of one orientation joined across gaps <= MAX_GAP_BP, gap bases adding no support."""
    groups = []
    for start, end, _exact in sorted(spans):
        # a wider gap starts a new group
        if not groups or start - groups[-1]['end'] > MAX_GAP_BP:
            groups.append(
                dict(
                    start=start,
                    end=end,
                    orientation=orientation,
                    supported_bp=end - start,
                    bridged_gaps=[],
                )
            )
            continue

        group = groups[-1]
        if start > group['end']:
            group['bridged_gaps'].append([group['end'], start])
            group['supported_bp'] += end - start
        elif end > group['end']:
            group['supported_bp'] += end - group['end']
        group['end'] = max(group['end'], end)

    for group in groups:
        group.update(
            support_start=group['start'],
            support_end=group['end'],
            coherent=True,
            bridged_bp=sum(gap_end - gap_start for gap_start, gap_end in group['bridged_gaps']),
        )
    return groups


def find_arrays(seq, motif):
    """Repeat-covered mask, arrays (>= MIN_UNITS units of support) and shorter components."""
    covered = bytearray(len(seq))
    arrays = []
    short = []
    for orientation, unit in (('F', motif), ('R', rc(motif))):
        spans = list(matching_spans(seq, rotations(unit)))
        for start, end, _exact in spans:
            covered[start:end] = b'\1' * (end - start)
        for group in components(spans, orientation):
            if group['supported_bp'] >= MIN_UNITS * len(motif):
                arrays.append(group)
            else:
                short.append(group)
    arrays.sort(key=lambda array: (array['start'], array['end'], array['orientation']))
    return covered, arrays, short


def forgive_edge_fragments(residual, arrays, short, unit_bp):
    """Residual repeat mask without array-edge overhangs and single terminal fragments."""
    residual = bytearray(residual)

    # a <1-unit overhang of an overlapping short match at an array edge is boundary ambiguity
    for start, end in list(runs(residual)):
        if end - start >= unit_bp:
            continue
        for array in arrays:
            # the run lies within unit_bp - 1 bp outside this array's start or end
            at_start = array['start'] - unit_bp + 1 <= start < end <= array['start']
            at_end = array['end'] <= start < end <= array['end'] + unit_bp - 1

            # and inside a short component that overlaps this array
            overlapping_short = any(
                component['start'] <= start
                and component['end'] >= end
                and component['start'] < array['end']
                and component['end'] > array['start']
                for component in short
            )
            if (at_start or at_end) and overlapping_short:
                residual[start:end] = b'\0' * (end - start)
                break

    # forgive a lone <=1-unit fragment at each intron edge, within EDGE_BP and outside the arrays
    if arrays:
        first = min(array['start'] for array in arrays)
        last = max(array['end'] for array in arrays)
        for side in ('donor', 'acceptor'):
            if side == 'donor':
                fragments = [run for run in runs(residual) if run[1] <= min(EDGE_BP, first)]
            else:
                fragments = [
                    run for run in runs(residual) if run[0] >= max(len(residual) - EDGE_BP, last)
                ]
            if len(fragments) == 1:
                start, end = fragments[0]
                if end - start <= unit_bp:
                    residual[start:end] = b'\0' * (end - start)
    return residual


def decide_architecture(seq, unit_bp, covered, arrays, unexplained_bp):
    """Architecture, reason and linker of the arrays. Any problem leaves mixed_unclear."""
    decision = dict(
        architecture='mixed_unclear',
        reason='',
        linker_start=None,
        linker_end=None,
        longest_nonrepeat_linker=0,
    )

    problems = []
    if set(seq) - set('ACGT'):
        problems.append('ambiguous_sequence')
    if unexplained_bp:
        problems.append('repeat_islands_outside_clean_arrays')
    if not arrays:
        problems.append('no_qualifying_array')
    if len(arrays) > 2:
        problems.append('more_than_two_arrays')
    if problems:
        decision['reason'] = ';'.join(sorted(problems))
        return decision

    if len(arrays) == 1:
        decision['architecture'] = arrays[0]['orientation']
        return decision

    # two arrays from here on, which must face opposite ways
    left, right = arrays
    if left['orientation'] == right['orientation']:
        decision['reason'] = 'separated_arrays_of_same_orientation'
        return decision
    joined = left['orientation'] + '-' + right['orientation']

    # an overlap (negative gap) must be under one unit, right array not nested
    gap = right['start'] - left['end']
    if gap < 0:
        if -gap >= unit_bp or right['end'] <= left['end']:
            decision['reason'] = 'ambiguous_opposite_array_overlap'
        else:
            decision['architecture'] = joined
        return decision

    # longest run of non-repeat ACGT bases between the arrays
    longest = run = 0
    for i in range(left['end'], right['start']):
        run = run + 1 if not covered[i] and seq[i] in 'ACGT' else 0
        longest = max(longest, run)
    decision['longest_nonrepeat_linker'] = longest

    # a linker >= MIN_LINKER_BP gives X-L-Y, a gap <= MAX_GAP_BP gives X-Y
    if longest >= MIN_LINKER_BP:
        decision.update(
            architecture=left['orientation'] + '-L-' + right['orientation'],
            linker_start=left['end'],
            linker_end=right['start'],
        )
    elif gap <= MAX_GAP_BP:
        decision['architecture'] = joined
    else:
        decision['reason'] = 'interarray_gap_without_qualifying_linker'
    return decision


def classify(seq, motif):
    """Tier, repeat coverage, arrays and architecture of one transcript-oriented intron."""
    seq = seq.upper()
    motif = motif.upper()
    if not seq or not motif or set(motif) - set('ACGT'):
        raise ValueError('Require a nonempty sequence and resolved motif')
    covered, arrays, short = find_arrays(seq, motif)

    # residual: repeat-covered bases outside every array, before and after edge forgiveness
    explained = bytearray(len(seq))
    for array in arrays:
        explained[array['start'] : array['end']] = b'\1' * (array['end'] - array['start'])
    residual = bytearray(cov and not exp for cov, exp in zip(covered, explained))
    unexplained = forgive_edge_fragments(residual, arrays, short, len(motif))

    covered_bp = sum(covered)
    result = dict(
        tier=tier(covered_bp, len(seq)),
        covered_bp=covered_bp,
        repeat_fraction=covered_bp / len(seq),
        arrays=arrays,
        repeat_bp_outside_arrays=sum(residual),
        unexplained_repeat_bp=sum(unexplained),
        ignored_edge_repeat_bp=sum(residual) - sum(unexplained),
        bridged_gap_bp=sum(array['bridged_bp'] for array in arrays),
    )
    result.update(decide_architecture(seq, len(motif), covered, arrays, sum(unexplained)))
    return result
