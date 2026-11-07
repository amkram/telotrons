"""Putative target-site duplications (TSDs) at primary telotrons.

Re-measures each primary locus from the genome, checks it against the catalogue and writes the
per-locus calls and totals per lineage and genome.
"""

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'lib'))
import measure_junctions as junctions
from _common import find_genome_files, open_maybe_gz, read_tsv, write_tsv

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--catalogue', required=True, help='stage 1 catalogue; tier=standard is tested')
parser.add_argument('--config', default=str(REPO / 'config/genomes.tsv'))
parser.add_argument('--refseq-dir', default=str(REPO / 'data/raw/refseq'))
parser.add_argument('--tara-dir', default=str(REPO / 'data/raw/tara'))
parser.add_argument('--out', default=str(REPO / 'work/03_tsd'))
args = parser.parse_args()

OUT = Path(args.out)
OUT.mkdir(parents=True, exist_ok=True)


def read_contigs(path, wanted):
    """The wanted contigs of a FASTA, reading only until all are found."""
    contigs = {}
    seqid = None
    chunks = []
    with open_maybe_gz(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if seqid in wanted:
                    contigs[seqid] = ''.join(chunks).upper()
                seqid = line[1:].split()[0]
                chunks = []
                if len(contigs) == len(wanted):
                    break
            elif seqid in wanted:
                chunks.append(line.strip())
        else:
            if seqid in wanted:
                contigs[seqid] = ''.join(chunks).upper()
    assert wanted <= contigs.keys(), ('missing_contigs', wanted - contigs.keys())
    return contigs


def locus_sequences(row, contigs, gaps):
    """Intron and its exon flanks (up to 60 bp, within the bounding exons) in transcript order."""
    strand = row['strand'] if row['strand'] in ('+', '-') else row['display_strand']
    key = (row['seqid'], int(row['start']), int(row['end']), strand)
    support = gaps.get(key, [])
    # a curated strand flip uses the opposite-strand models with exon sides swapped
    if not support and row['tsd_annotation_support_basis'] == 'established_strand_correction':
        opposite = gaps.get(key[:3] + ('-' if strand == '+' else '+',), [])
        support = [(model, gene, right, left, kind) for model, gene, left, right, kind in opposite]

    flanked = junctions.window(contigs, key)
    assert flanked is not None, (row['locus_id'], 'missing_window')
    left, intron, right = junctions.bound_flanks(flanked, support)
    assert bool(support) == bool(int(row['tsd_annotation_supported'])), row['locus_id']
    assert hashlib.sha256(intron.encode()).hexdigest() == row['intron_sha256'], row['locus_id']
    return dict(
        locus_id=row['locus_id'],
        genome_id=row['genome_id'],
        seqid=key[0],
        start=key[1],
        end=key[2],
        strand=key[3],
        motif=row['motif'],
        left=left,
        intron=intron,
        right=right,
        annotation_supported=bool(support),
    )


lineage = {row['genome_id']: row['lineage'] for row in read_tsv(args.config)}
primary = [row for row in read_tsv(args.catalogue) if row['tier'] == 'standard']
assert primary and len({row['locus_id'] for row in primary}) == len(primary)

loci = []
records = []
for genome_id in sorted({row['genome_id'] for row in primary}):
    rows = [row for row in primary if row['genome_id'] == genome_id]
    fasta, gff = find_genome_files(genome_id, args.refseq_dir, args.tara_dir)
    assert fasta and gff, (genome_id, 'missing_genome_files')
    gaps = junctions.annotation_gaps(gff)
    contigs = read_contigs(fasta, {row['seqid'] for row in rows})

    for row in rows:
        record = locus_sequences(row, contigs, gaps)
        records.append(record)
        observed = junctions.measure(
            record['left'], record['intron'], record['right'], record['motif']
        )
        for field in ('core_score', 'core_bp', 'core_back', 'core_fwd', 'core_sequence'):
            assert str(observed[field]) == row['tsd_' + field], (row['locus_id'], field)
        putative_tsd = int(record['annotation_supported'] and observed['core_score'] > 0)
        assert putative_tsd == int(row['tsd_adjustment_used']), row['locus_id']
        loci.append(
            dict(
                locus_id=row['locus_id'],
                genome_id=genome_id,
                lineage=lineage[genome_id],
                motif=row['motif'],
                core_score=int(observed['core_score']),
                core_bp=int(observed['core_bp']),
                putative_TSD=putative_tsd,
                left_bp=len(record['left']),
                right_bp=len(record['right']),
            )
        )
    print(genome_id, len(loci), 'loci', flush=True)

write_tsv(OUT / 'TSD_loci.tsv', loci)
with gzip.open(OUT / 'locus_sequences.jsonl.gz', 'wt') as handle:
    for record in records:
        handle.write(json.dumps(record, separators=(',', ':')) + '\n')

# TSD counts over all loci, per lineage and per genome
groups = (
    [('all', 'ALL')]
    + [('lineage', name) for name in sorted({row['lineage'] for row in loci})]
    + [('genome_id', name) for name in sorted({row['genome_id'] for row in loci})]
)
summary = []
for group_type, group in groups:
    members = [row for row in loci if group_type == 'all' or row[group_type] == group]
    tsds = sum(row['putative_TSD'] for row in members)
    summary.append(
        dict(
            group_type=group_type,
            group=group,
            loci=len(members),
            TSD=tsds,
            percent=100 * tsds / len(members),
        )
    )
write_tsv(OUT / 'TSD_summary.tsv', summary)
print(json.dumps(summary[0]))
