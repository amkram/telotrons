"""Architecture-specific coverage: strict >85% except linker classes >50%."""

NONLINKER = frozenset(['F', 'R', 'F-R', 'R-F'])
LINKER = frozenset(['F-L-R', 'R-L-F'])


def apply_coverage_gate(architecture, covered_bp, intron_bp):
    """Architecture after the coverage threshold of its class, with the threshold and outcome."""
    assert 0 <= covered_bp <= intron_bp and intron_bp > 0
    assert architecture in NONLINKER | LINKER | {'mixed_unclear'}

    # mixed_unclear has no threshold, so it always passes
    if architecture in NONLINKER:
        threshold = 85
    elif architecture in LINKER:
        threshold = 50
    else:
        threshold = None
    passed = threshold is None or covered_bp * 100 > threshold * intron_bp
    reason = '' if passed else f'{architecture}_requires_strictly_gt{threshold}_percent_coverage'
    return dict(
        architecture=architecture if passed else 'mixed_unclear',
        architecture_coverage_threshold=threshold,
        architecture_coverage_pass=int(passed),
        coverage_gate_reason=reason,
    )
