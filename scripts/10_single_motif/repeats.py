"""Motif families of period 1-8 and a ctypes wrapper around the C repeat scanner (repeats.c)."""

import ctypes
import itertools

BASES = 'ACGT'


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def units(motif):
    """Every rotation of the motif and of its reverse complement."""
    strands = [motif, reverse_complement(motif)]
    return sorted({strand[i:] + strand[:i] for strand in strands for i in range(len(motif))})


def primitive(word):
    """True unless the word repeats a shorter word."""
    return not any(
        len(word) % k == 0 and word == word[:k] * (len(word) // k) for k in range(1, len(word))
    )


def canonical(motif):
    return min(units(motif))


class Hit(ctypes.Structure):
    _fields_ = [('id', ctypes.c_int32), ('screen', ctypes.c_int32), ('tract', ctypes.c_int32)]


class Scanner:
    """Screen and tract coverage of every primitive motif family, computed in C."""

    def __init__(self, library):
        # canonical motif id of every word of length 1 to 8 (-1 if not primitive)
        self.motifs = []
        self.ids = {}
        word_ids = []
        offsets = [0] * 9
        for k in range(1, 9):
            offsets[k] = len(word_ids)
            for letters in itertools.product(BASES, repeat=k):
                word = ''.join(letters)
                if not primitive(word):
                    word_ids.append(-1)
                    continue
                motif = canonical(word)
                if motif not in self.ids:
                    self.ids[motif] = len(self.ids)
                    self.motifs.append(motif)
                word_ids.append(self.ids[motif])
        self.patterns = [''.join(units(motif)).encode() for motif in self.motifs]

        # C function signatures
        self.lib = ctypes.CDLL(str(library))
        self.lib.create_scanner.argtypes = [
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32),
            ctypes.c_int32,
        ]
        self.lib.create_scanner.restype = ctypes.c_void_p
        self.lib.scan.argtypes = [
            ctypes.c_void_p,
            ctypes.c_char_p,
            ctypes.c_int32,
            ctypes.POINTER(Hit),
            ctypes.c_int,
        ]
        self.lib.scan.restype = ctypes.c_int
        self.lib.target_coverage.argtypes = [
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_int),
        ]
        self.lib.target_tract.argtypes = self.lib.target_coverage.argtypes
        self.lib.free_scanner.argtypes = [ctypes.c_void_p]

        self.handle = self.lib.create_scanner(
            len(self.ids),
            (ctypes.c_int32 * len(word_ids))(*word_ids),
            (ctypes.c_int32 * 9)(*offsets),
            len(word_ids),
        )
        if not self.handle:
            raise MemoryError('repeat scanner')
        # reusable output buffers
        self.hits = (Hit * len(self.ids))()
        self.out = (ctypes.c_int * 3)()

    def scan(self, sequence, threshold=50):
        """(motif id, screen bp, tract bp) of every motif where either exceeds threshold %."""
        if isinstance(sequence, str):
            sequence = sequence.encode()
        n_hits = self.lib.scan(self.handle, sequence, len(sequence), self.hits, threshold)
        return [(hit.id, hit.screen, hit.tract) for hit in self.hits[:n_hits]]

    def coverage(self, sequence, motif_id):
        """(screen-rule bp, exact-unit bp with gaps of <= 3 bp closed) of one motif."""
        if isinstance(sequence, str):
            sequence = sequence.encode()
        period = len(self.motifs[motif_id])
        pattern = self.patterns[motif_id]
        self.lib.target_coverage(
            sequence, len(sequence), pattern, period, len(pattern) // period, self.out
        )
        return tuple(self.out[:2])

    def tract(self, sequence, motif_id):
        """(tract bp, tract bp in the first 12 bp, tract bp in the last 12 bp) of one motif."""
        if isinstance(sequence, str):
            sequence = sequence.encode()
        period = len(self.motifs[motif_id])
        pattern = self.patterns[motif_id]
        self.lib.target_tract(
            sequence, len(sequence), pattern, period, len(pattern) // period, self.out
        )
        return tuple(self.out)

    def close(self):
        if self.handle:
            self.lib.free_scanner(self.handle)
            self.handle = None
