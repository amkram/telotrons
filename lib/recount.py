"""Recount one genome's raw-scan ITS components outside the projected gene spans as well."""

import collections, json
from pathlib import Path
from audit_its_opportunity import (
    count_acgt,
    fasta,
    gene_mask,
    merge,
    overlaps,
    read_tsv,
    subtract,
    write_tsv,
)
from run_projection import alignments

# a projected span counts with >= 50% query coverage, >= 40% identity, >= 100 aligned residues
# and a score >= 80% of the query's best
MIN_COVERAGE = 0.5
MIN_IDENTITY = 0.4
MIN_ALIGNED = 100
MIN_RELATIVE_SCORE = 0.8


def projected_mask(paths):
    """Merged spans of every alignment passing the filters, per contig."""
    # best score per query over all alignments
    best = {}
    for path in paths:
        for alignment in alignments(path):
            query = alignment['query']
            best[query] = max(best.get(query, -(10**30)), alignment['score'])

    spans = collections.defaultdict(list)
    for path in paths:
        for alignment in alignments(path):
            aligned = alignment['qend'] - alignment['qstart']
            top = best[alignment['query']]
            if (
                top > 0
                and alignment['score'] > 0
                and aligned / alignment['qlen'] >= MIN_COVERAGE
                and alignment['identity'] >= MIN_IDENTITY
                and aligned >= MIN_ALIGNED
                and alignment['score'] >= MIN_RELATIVE_SCORE * top
            ):
                spans[alignment['seqid']].append((alignment['start0'], alignment['end0']))
    return {contig: merge(pieces) for contig, pieces in spans.items()}


def masked_acgt(sequence, intervals):
    """A/C/G/T bases inside disjoint intervals of a sequence."""
    return sum(count_acgt(sequence[start:end]) for start, end in intervals)


def eligible_acgt(fasta_path, annotation, caps, projected):
    """Eligible A/C/G/T bp outside the gene mask and caps, without and with projected spans."""
    before = 0
    after = 0
    for contig, sequence in fasta(fasta_path):
        # annotation plus cap mask, merged, then projected bp not already in it
        base = merge(annotation.get(contig, []) + caps.get(contig, []))
        extra = subtract(projected.get(contig, []), base)
        remaining = count_acgt(sequence) - masked_acgt(sequence, base)
        before += remaining
        after += remaining - masked_acgt(sequence, extra)
    return before, after


def recount(genome, fasta_path, gff_path, raw_folder, folder, paths, known):
    """Count raw-scan components outside the projected spans into folder/summary.tsv."""
    genome_id = genome['genome_id']
    raw = Path(raw_folder) / genome_id
    scan = json.loads(Path(f'{raw}.json').read_text())
    kept = [
        (row['seqid'], int(row['start0']), int(row['end0']))
        for row in read_tsv(f'{raw}.components.tsv.gz')
    ]
    assert len(kept) == scan['n_components']

    # gene, cap and projected masks. Without projection the eligible bp must match the raw scan
    caps = collections.defaultdict(list)
    for row in read_tsv(f'{raw}.caps.tsv.gz'):
        caps[row['seqid']].append((int(row['start0']), int(row['end0'])))
    projected = projected_mask(paths)
    before, eligible = eligible_acgt(fasta_path, gene_mask(gff_path), caps, projected)
    assert before == scan['eligible_acgt_bp'], (genome_id, before, scan['eligible_acgt_bp'])

    # no counted component may overlap a known telotron (0-based, controls have none)
    telotrons = collections.defaultdict(list)
    for contig, start, end in known:
        telotrons[contig].append((start, end))
    telotrons = {contig: merge(pieces) for contig, pieces in telotrons.items()}
    hits = sum(overlaps(start, end, telotrons.get(contig, [])) for contig, start, end in kept)
    assert hits == 0, (genome_id, 'known telotrons overlap counted ITS components', hits)

    # the reported count: components outside the projected spans too
    n_components = sum(
        not overlaps(start, end, projected.get(contig, [])) for contig, start, end in kept
    )
    summary = dict(
        genome_id=genome_id,
        role=genome['role'],
        clade=genome['lineage'],
        assembly_class=genome['assembly_class'],
        organism=genome['label'],
        n_components=n_components,
        eligible_acgt_bp=eligible,
    )
    write_tsv(Path(folder) / 'summary.tsv', [summary])
    print(genome_id, n_components, 'components;', eligible, 'eligible A/C/G/T bp', flush=True)
