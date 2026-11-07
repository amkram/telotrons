"""TSD-aware 8-bp endpoint gate for clean array shapes (0-based half-open, transcript-oriented)."""

EDGE_BP = 8


def boundary_reach(has_arrays, donor_ok, acceptor_ok):
    """Which intron ends the arrays reach: none, both, donor_only, acceptor_only or neither."""
    if not has_arrays:
        return 'none'
    if donor_ok and acceptor_ok:
        return 'both'
    if donor_ok:
        return 'donor_only'
    if acceptor_ok:
        return 'acceptor_only'
    return 'neither'


def apply_boundary_gate(shape, intron_bp, arrays, core, annotation_supported):
    """Shape kept only when its outermost arrays lie within EDGE_BP of the expected intron edges."""
    # TSD = measured repeat-excluded copy core, not an ancestral TSD length
    use_tsd = bool(annotation_supported and core['core_score'] > 0)

    # expected array span, moved inward past the copy core when used
    start = int(core['core_fwd']) if use_tsd else 0
    end = intron_bp - (int(core['core_back']) if use_tsd else 0)
    assert 0 <= start <= end <= intron_bp

    # outermost repeat support across all arrays
    first = min((array['support_start'] for array in arrays), default=None)
    last = max((array['support_end'] for array in arrays), default=None)

    # distances to the expected edges, raw_* tested against the unshifted intron ends
    left_distance = None if first is None else abs(first - start)
    right_distance = None if last is None else abs(end - last)
    left_ok = left_distance is not None and left_distance <= EDGE_BP
    right_ok = right_distance is not None and right_distance <= EDGE_BP
    raw_left_ok = first is not None and first <= EDGE_BP
    raw_right_ok = last is not None and intron_bp - last <= EDGE_BP

    problems = []
    if not left_ok:
        problems.append('donor_array_not_within_8bp_of_intron_or_tsd_boundary')
    if not right_ok:
        problems.append('acceptor_array_not_within_8bp_of_intron_or_tsd_boundary')
    return dict(
        shape_architecture=shape,
        architecture=shape if left_ok and right_ok else 'mixed_unclear',
        boundary_gate_reason=';'.join(problems),
        edge_tolerance_bp=EDGE_BP,
        tsd_adjustment_used=int(use_tsd),
        tsd_boundary_basis='measured_repeat_excluded_copy_core' if use_tsd else 'intron_boundary',
        tsd_annotation_supported=int(bool(annotation_supported)),
        donor_boundary_distance_bp=left_distance,
        acceptor_boundary_distance_bp=right_distance,
        donor_boundary_pass=int(left_ok),
        acceptor_boundary_pass=int(right_ok),
        boundary_reach=boundary_reach(bool(arrays), left_ok, right_ok),
        # clean only because the TSD shift moved the expected edges
        clean_depends_on_tsd_adjustment=int(
            shape != 'mixed_unclear' and left_ok and right_ok and not (raw_left_ok and raw_right_ok)
        ),
    )
