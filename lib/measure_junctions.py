"""Repeat-excluded junction copy cores (present-day exact copies, not ancestral TSDs)."""

from collections import defaultdict
import gzip
import math


from _common import rc
from repeat_policy import coverage, motif_units

CAP = 60
MIN_CORE = 6
MIN_NONSPLICE = 4


def informative(seq, back, fwd):
    """Whether a copy core is long and complex enough to count."""
    # long enough, with enough copy bases beyond the two splice-site bases per side
    if len(seq) < MIN_CORE or max(0, back - 2) + max(0, fwd - 2) < MIN_NONSPLICE:
        return False

    # low-complexity filters: ambiguity, base diversity, composition, entropy
    if set(seq) - set('ACGT') or len(set(seq)) < 3:
        return False
    freqs = [seq.count(base) / len(seq) for base in 'ACGT' if base in seq]
    if max(freqs) > 0.75 or -sum(freq * math.log2(freq) for freq in freqs) < 1.2:
        return False

    # not a tandem repeat of period 1, 2 or 3
    return not any(seq == (seq[:period] * len(seq))[: len(seq)] for period in (1, 2, 3))


def side_span(host, intronic, host_mask, intron_mask, reverse=False, cap=CAP):
    """Length of the exact host/intron copy outward from the junction, and why it stopped."""
    # reverse=True measures leftwards by reversing every input
    if reverse:
        host = host[::-1]
        intronic = intronic[::-1]
        host_mask = host_mask[::-1]
        intron_mask = intron_mask[::-1]

    limit = min(cap, len(host), len(intronic))
    for i in range(limit):
        if host[i] not in 'ACGT' or intronic[i] not in 'ACGT':
            return i, 'ambiguous_base'
        if host[i] != intronic[i]:
            return i, 'mismatch'
        if host_mask[i]:
            return i, 'repeat_in_exonic_copy'
        if intron_mask[i]:
            return i, 'repeat_in_intronic_copy'
    return limit, ('search_cap' if limit == cap else 'sequence_boundary')


def copy_core(left, intron, right, masks, cap=CAP):
    """Copy core spanning both junctions, with flank and intron repeat bases excluded."""
    left_mask, intron_mask, right_mask = masks

    # back: left flank end vs intron end, fwd: right flank start vs intron start
    back, back_stop = side_span(left, intron, left_mask, intron_mask, True, cap)
    fwd, fwd_stop = side_span(right, intron, right_mask, intron_mask, False, cap)
    sequence = (left[-back:] if back else '') + intron[:fwd]
    assert sequence == (intron[-back:] if back else '') + right[:fwd]

    # a copy longer than the intron overlaps itself and scores 0
    overlap = back + fwd > len(intron)
    return dict(
        core_back=back,
        core_fwd=fwd,
        core_bp=back + fwd,
        core_sequence=sequence,
        core_back_stop=back_stop,
        core_fwd_stop=fwd_stop,
        core_overlap=overlap,
        core_score=back + fwd if not overlap and informative(sequence, back, fwd) else 0,
        core_nonsplice_bp=max(0, back - 2) + max(0, fwd - 2),
        # either side stopped at the search cap or the sequence end
        core_censored=back_stop in ('search_cap', 'sequence_boundary')
        or fwd_stop in ('search_cap', 'sequence_boundary'),
    )


def measure(left, intron, right, motif, cap=CAP):
    """Copy core of one intron under the strict motif repeat masks."""
    units = motif_units(motif)
    masks = [coverage(seq, units) for seq in (left, intron, right)]
    return copy_core(left, intron, right, masks, cap)


def annotation_gaps(path):
    """1-based closed gaps between consecutive blocks of each parent: exons, else CDS."""
    blocks = defaultdict(lambda: defaultdict(set))
    parent_of = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            # skip header lines, stopping at an embedded FASTA section
            if line.startswith('#'):
                if line.strip() == '##FASTA':
                    break
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9:
                continue

            # map each ID to its first Parent (transcript -> gene)
            attributes = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            if 'ID' in attributes and 'Parent' in attributes:
                parent_of[attributes['ID']] = attributes['Parent'].split(',')[0]

            # exon/CDS coordinates per parent transcript
            if fields[2] not in ('exon', 'CDS'):
                continue
            for parent in attributes.get('Parent', '').split(','):
                if parent:
                    blocks[parent][fields[2]].add(
                        (fields[0], fields[6], int(fields[3]), int(fields[4]))
                    )

    # gaps between consecutive non-abutting blocks, keyed by (contig, start, end, strand)
    gaps = defaultdict(list)
    for parent, kinds in blocks.items():
        kind = 'exon' if kinds.get('exon') else 'CDS'
        spans_by_strand = defaultdict(list)
        for seqid, strand, start, end in kinds[kind]:
            spans_by_strand[seqid, strand].append((start, end))
        for (seqid, strand), spans in spans_by_strand.items():
            spans.sort()
            for upstream, downstream in zip(spans, spans[1:]):
                if downstream[0] <= upstream[1] + 1:
                    continue

                # flanking block lengths, upstream first in transcript orientation
                left_bp = upstream[1] - upstream[0] + 1
                right_bp = downstream[1] - downstream[0] + 1
                if strand == '-':
                    left_bp, right_bp = right_bp, left_bp
                gaps[seqid, upstream[1] + 1, downstream[0] - 1, strand].append(
                    (parent, parent_of.get(parent, parent), left_bp, right_bp, kind)
                )
    return gaps


def window(contigs, key, flank=CAP):
    """Left flank, intron and right flank of a 1-based closed intron, in transcript orientation."""
    seqid, start, end, strand = key
    chrom = contigs.get(seqid, '')
    if not 1 <= start <= end <= len(chrom):
        return None
    left = chrom[max(0, start - 1 - flank) : start - 1]
    intron = chrom[start - 1 : end]
    right = chrom[end : end + flank]

    # minus strand: reverse-complement and swap flanks
    if strand == '-':
        left, intron, right = rc(right), rc(intron), rc(left)
    return left, intron, right


def bound_flanks(flanked, support):
    """Trim the flanks to the shortest supporting exon on each side."""
    if not support:
        return flanked
    left_bp = min(model[2] for model in support)
    right_bp = min(model[3] for model in support)
    return flanked[0][-left_bp:], flanked[1], flanked[2][:right_bp]
