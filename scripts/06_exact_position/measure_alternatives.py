"""Repeat coverage and composition state of each alternative intron under the focal motif."""

import hashlib
import json

from exact_position import OUT
from paired_repeat_tracts import tracts
from repeat_policy import coverage, motif_units

TELOTRON_MIN = 0.5  # repeat fraction above this is TELOTRON
ORDINARY_MAX = 0.3  # below this, with no paired tract of ORDINARY_MAX_UNITS units, ORDINARY
ORDINARY_MAX_UNITS = 4


def composition_state(fraction, longest, motif):
    """TELOTRON, ORDINARY or MIXED_OR_AMBIGUOUS from coverage and the longest paired tract."""
    if fraction > TELOTRON_MIN:
        return 'TELOTRON'
    if fraction < ORDINARY_MAX and longest < ORDINARY_MAX_UNITS * len(motif):
        return 'ORDINARY'
    return 'MIXED_OR_AMBIGUOUS'


def main():
    rows = json.loads((OUT / 'alternative_sequences.json').read_text())
    for row in rows:
        seq = row['sequence']
        motif = row['motif']
        assert hashlib.sha256(seq.encode()).hexdigest() == row['sequence_sha256']
        fraction = sum(coverage(seq, motif_units(motif))) / len(seq)
        longest = max(
            (end - start for start, end, *_ in tracts(seq, motif, min_units=2)), default=0
        )
        row.update(
            repeat_fraction=fraction,
            longest_paired_tract=longest,
            state_by_focal_motif=composition_state(fraction, longest, motif),
            splice=seq[:2] + '-' + seq[-2:],
        )

    n_telomeric = sum(row['state_by_focal_motif'] == 'TELOTRON' for row in rows)
    (OUT / 'alternative_composition.json').write_text(json.dumps(rows, indent=2) + '\n')
    print(len(rows), 'alternatives,', n_telomeric, 'TELOTRON by the focal motif', flush=True)


if __name__ == '__main__':
    main()
