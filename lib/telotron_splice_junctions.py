"""Count RNA-seq alignments whose N skip matches a locus' intron boundaries exactly.

Counting policy for telotrons and ordinary introns alike: primary alignment records with
MAPQ >= 20 and at least 10 contiguous aligned bases on each side of the skip. Counts are alignment
records, not independent fragments or molecules.
"""

import csv
from rna_evidence import atomic_path, bam_header, checked_output, junctions, primary_records

MIN_MAPQ = 20
MIN_ANCHOR = 10
COLUMNS = ['seqid', 'start', 'end', 'reads_over_locus', 'exact_junction_reads']


def check_loci(loci, contigs):
    """1-based inclusive loci must fit their BAM reference and appear once."""
    seen = set()
    for locus in loci:
        chrom, start, end = locus['seqid'], int(locus['start']), int(locus['end'])
        if chrom not in contigs or not 1 <= start <= end <= contigs[chrom]:
            raise ValueError(f'Locus outside BAM references: {chrom}:{start}-{end}')
        key = chrom, start, end, locus.get('strand', '')
        if key in seen:
            raise ValueError(f'Duplicate locus: {key}')
        seen.add(key)


def count_locus(bam, seqid, start, end):
    """(reads over the locus, reads with an N skip exactly from start to end)."""
    # 0x904: unmapped, secondary, supplementary
    command = ['samtools', 'view', '-F', '0x904', '-q', str(MIN_MAPQ), str(bam)]
    records = checked_output([*command, f'{seqid}:{start}-{end}'])
    overlapping = exact = 0
    for fields in primary_records(records, MIN_MAPQ):
        overlapping += 1
        if 'N' not in fields[5]:
            continue
        # a read counts once, at its first skip matching both boundaries
        for skip_start, skip_end, left, right in junctions(int(fields[3]), fields[5]):
            if skip_start == start and skip_end == end and min(left, right) >= MIN_ANCHOR:
                exact += 1
                break
    return overlapping, exact


def count_junctions(bam, loci, output):
    """Write one row per locus (seqid, start, end, read counts) to output and return the rows."""
    check_loci(loci, bam_header(bam))
    rows = []
    for locus in loci:
        seqid, start, end = locus['seqid'], int(locus['start']), int(locus['end'])
        rows.append([seqid, start, end, *count_locus(bam, seqid, start, end)])

    with atomic_path(output) as temporary:
        with open(temporary, 'w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(COLUMNS)
            writer.writerows(rows)
    return rows
