"""Telomeric tracts from the C scanner, their checks, merged components and terminal caps."""

import ctypes, os, subprocess, sys
from pathlib import Path
from audit_its_opportunity import merge

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


class Tract(ctypes.Structure):
    """One tract as the C scanner returns it."""

    _fields_ = [
        ('start', ctypes.c_int64),
        ('end', ctypes.c_int64),
        ('pattern', ctypes.c_int32),
        ('mismatches', ctypes.c_int32),
        ('exact_units', ctypes.c_int32),
    ]


def revcomp(sequence):
    """Reverse complement."""
    return sequence.translate(COMPLEMENT)[::-1]


def build(folder):
    """Compile lib/its_scan_core.c into folder (when missing or older) and load it."""
    source = Path(__file__).resolve().parent / 'its_scan_core.c'
    library = Path(folder) / 'its_scan_core.so'
    library.parent.mkdir(parents=True, exist_ok=True)

    # compile to a private name, then rename, so parallel jobs never load a partial file
    if not library.exists() or source.stat().st_mtime > library.stat().st_mtime:
        partial = library.with_name(f'{library.name}.{os.getpid()}')
        shared = ['-dynamiclib'] if sys.platform == 'darwin' else ['-shared', '-fPIC']
        command = ['cc', '-O3', '-std=c11', '-Wall', '-Wextra', '-Werror', *shared]
        subprocess.run([*command, str(source), '-o', str(partial)], check=True)
        os.replace(partial, library)

    # ctypes signatures of the C functions
    core = ctypes.CDLL(str(library))
    core.scan_its.argtypes = [
        ctypes.c_char_p,
        ctypes.c_int64,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.POINTER(Tract)),
    ]
    core.scan_its.restype = ctypes.c_int64
    core.free_its.argtypes = [ctypes.POINTER(Tract)]
    return core


def patterns(motif):
    """Every rotation of the motif and of its reverse complement."""
    if not 2 <= len(motif) <= 8 or set(motif) - set('ACGT'):
        raise ValueError('Unsupported motif')
    strands = sorted([motif, revcomp(motif)])
    return [strand[i:] + strand[:i] for strand in strands for i in range(len(strand))]


def scan(sequence, motif, core):
    """Maximal phased tracts (start, end, pattern, mismatches, exact_units) of a sequence."""
    # the C side allocates the tract array, freed here after copying
    tracts = ctypes.POINTER(Tract)()
    n_tracts = core.scan_its(
        sequence.encode(),
        len(sequence),
        ''.join(patterns(motif)).encode(),
        len(motif),
        ctypes.byref(tracts),
    )
    if n_tracts < 0:
        raise RuntimeError(f'C scanner failed: {n_tracts}')
    try:
        fields = [name for name, _ in Tract._fields_]
        return [tuple(getattr(tract, name) for name in fields) for tract in tracts[:n_tracts]]
    finally:
        core.free_its(tracts)


def validate_tracts(sequence, motif, tracts):
    """Assert the C scanner's tract rules for every tract."""
    phases = patterns(motif)
    unit = len(motif)
    for start, end, pattern_index, mismatches, exact in tracts:
        # whole units, at least four, at least one of them exact
        assert 0 <= start < end <= len(sequence)
        assert (end - start) % unit == 0 and end - start >= 4 * unit and exact > 0
        pattern = phases[pattern_index]
        assert not (set(sequence[start:end]) - set('ACGT'))

        # mismatches of each unit against the tract's phase pattern
        distances = [
            sum(base != expected for base, expected in zip(sequence[i : i + unit], pattern))
            for i in range(start, end, unit)
        ]
        # at most one substitution per adjacent pair, totals equal to the C output
        assert all(left + right <= 1 for left, right in zip(distances, distances[1:]))
        assert max(distances) <= 1
        assert sum(distances) == mismatches
        assert distances.count(0) == exact

        # maximality: the unit just outside each end must not extend the tract
        for outside in [start - unit, end]:
            if 0 <= outside and outside + unit <= len(sequence):
                neighbour = sequence[outside : outside + unit]
                distance = sum(base != expected for base, expected in zip(neighbour, pattern))
                edge = distances[0] if outside < start else distances[-1]
                assert set(neighbour) - set('ACGT') or distance > 1 or distance + edge > 1


def components(tracts):
    """ITS components: overlapping or touching tracts merged (each tract spans >= 4 units)."""
    return merge((start, end) for start, end, *_ in tracts)


def cap_mask(sequence, motif, window=1000, fraction=0.15):
    """Terminal telomeric caps: end runs of windows whose motif bases reach the fraction."""
    length = len(sequence)
    unit = len(motif)
    reverse = revcomp(motif)

    # no self-overlap, so str.count matches tidk's occurrence counting
    assert not any(
        motif[:size] == motif[-size:] or reverse[:size] == reverse[-size:]
        for size in range(1, unit)
    )

    # per window: whether motif occurrences cover the fraction of its true width
    flags = []
    for start in range(0, length, window):
        part = sequence[start : start + window]
        count = part.count(motif) + part.count(reverse)
        flags.append(count * unit / len(part) >= fraction)

    # grow each end span while consecutive windows are flagged
    left = 0
    right = length
    for i, flag in enumerate(flags):
        if not flag:
            break
        left = min((i + 1) * window, length)
    for i in range(len(flags) - 1, -1, -1):
        if not flags[i]:
            break
        right = i * window

    # [span] * bool keeps an end span only when it is non-empty
    return merge([(0, left)] * (left > 0) + [(right, length)] * (right < length))
