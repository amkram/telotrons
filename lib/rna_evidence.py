"""Checked subprocess calls, BAM checks, CIGAR parsing and atomic writes for RNA evidence."""

from contextlib import contextmanager
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


def checked_output(command):
    """stdout of a command. A nonzero exit raises with its stderr."""
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {shlex.join(command)}\n"
            f"{result.stderr.strip()}"
        )
    return result.stdout


def parse_cigar(cigar):
    """(length, operation) pairs of a mapped read's CIGAR string."""
    if not re.fullmatch(r'(?:[1-9][0-9]*[MIDNSHP=X])+', cigar):
        raise ValueError(f'Malformed mapped-read CIGAR: {cigar!r}')
    return [(int(length), op) for length, op in re.findall(r'([0-9]+)([MIDNSHP=X])', cigar)]


def junctions(position, cigar):
    """(start, end, left, right) per N skip: 1-based inclusive skip, M/=/X anchor bp each side."""
    operations = parse_cigar(cigar)
    # reference position where the current operation starts
    reference = position
    for i, (length, operation) in enumerate(operations):
        if operation == 'N':
            # anchors: contiguous match runs just before and just after the skip
            left = right = 0
            for run, code in reversed(operations[:i]):
                if code not in 'M=X':
                    break
                left += run
            for run, code in operations[i + 1 :]:
                if code not in 'M=X':
                    break
                right += run
            yield reference, reference + length - 1, left, right

        # operations that consume reference bases
        if operation in 'MDN=X':
            reference += length


def primary_records(text, min_mapq=0, exclude_flags=0x904):
    """SAM records of mapped primary alignments with MAPQ >= min_mapq."""
    for line in text.splitlines():
        if not line or line.startswith('@'):
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            raise ValueError('Truncated SAM record')
        # default 0x904 drops unmapped, secondary and supplementary records
        if int(fields[1]) & exclude_flags or int(fields[4]) < min_mapq:
            continue
        # raises on a malformed CIGAR
        parse_cigar(fields[5])
        yield fields


def bam_header(bam):
    """Reference lengths of a BAM, after checking the file and its index."""
    checked_output(['samtools', 'quickcheck', '-v', str(bam)])
    header = checked_output(['samtools', 'view', '-H', str(bam)])

    # reference lengths from the @SQ header lines
    contigs = {}
    for line in header.splitlines():
        if line.startswith('@SQ\t'):
            fields = dict(item.split(':', 1) for item in line.split('\t')[1:] if ':' in item)
            contigs[fields['SN']] = int(fields['LN'])
    if not contigs:
        raise ValueError(f'No reference contigs in {bam}')

    # region query fails fast on an unindexed BAM
    chrom = next(iter(contigs))
    checked_output(['samtools', 'view', '-c', str(bam), f'{chrom}:1-1'])
    return contigs


@contextmanager
def atomic_path(output):
    """A temporary path beside output, renamed into place only on success."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f'.{output.name}.', dir=output.parent)
    os.close(descriptor)
    try:
        yield temporary
        os.replace(temporary, output)
    finally:
        # a failed write leaves no partial file behind
        if os.path.exists(temporary):
            os.unlink(temporary)
