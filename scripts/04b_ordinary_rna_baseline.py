#!/usr/bin/env python3
"""Exact junction reads at a seeded sample of 2,000 ordinary annotated introns of one genome of
config/rna_runs.tsv, with the same BAM and counting policy as its telotrons (04_rna_support.py).

Ordinary introns: annotated introns (lib/derive_introns.py, the set stage 1 screens) of >= 30 bp
with < 30% coverage by the lineage telomere motif (the stage 1 coverage measure), excluding
catalogue loci. Writes work/04_rna/<genome>/ordinary_sample.tsv, ordinary_junctions.tsv and
ordinary_summary.json.
"""

import argparse, json, random, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'lib'))
from _common import find_genome_files, rc
from audit_its_opportunity import fasta, read_tsv, rows_by_genome, write_tsv
from derive_introns import derive
from repeat_policy import coverage, motif_units
from telotron_splice_junctions import count_junctions

SEED = 20261105
SAMPLE_SIZE = 2000
MIN_LENGTH = 30
MAX_COVERAGE = 0.30

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--genome-id', required=True)
parser.add_argument('--catalogue', default=str(REPO / 'work/catalogue/catalogue.tsv.gz'))
parser.add_argument('--genomes', default=str(REPO / 'config/genomes.tsv'))
parser.add_argument('--runs', default=str(REPO / 'config/rna_runs.tsv'))
parser.add_argument('--refseq-dir', default=str(REPO / 'data/raw/refseq'))
parser.add_argument('--tara-dir', default=str(REPO / 'data/raw/tara'))
parser.add_argument('--out', default=str(REPO / 'work/04_rna'))
parser.add_argument('--bam', help='default: work/04_rna/<genome>/<genome>.splice_sr.bam')
args = parser.parse_args()

genome_id = args.genome_id
assert genome_id in rows_by_genome(args.runs), (genome_id, 'no RNA runs')
units = tuple(sorted(motif_units(rows_by_genome(args.genomes)[genome_id]['motif'])))
folder = Path(args.out) / genome_id
bam = args.bam or folder / f'{genome_id}.splice_sr.bam'

# annotated introns and the contigs that carry them
fasta_path, gff_path = find_genome_files(genome_id, args.refseq_dir, args.tara_dir)
assert fasta_path and gff_path, (genome_id, 'missing assembly or annotation under data/raw')
introns = derive(gff_path)
wanted = {contig for contig, start, end, strand in introns}
contigs = {contig: sequence for contig, sequence in fasta(fasta_path) if contig in wanted}

# this genome's telotron loci (0-based half-open), kept out of the ordinary pool
telotrons = {
    (row['seqid'], int(row['start']) - 1, int(row['end']))
    for row in read_tsv(args.catalogue)
    if row['genome_id'] == genome_id
}
assert telotrons, genome_id

# ordinary introns, each coordinate once, in coordinate order
pool = {}
for contig, start, end, strand in sorted(introns):
    key = (contig, start, end)
    if end - start < MIN_LENGTH or key in telotrons or key in pool:
        continue
    intron = contigs[contig][start:end]
    if strand == '-':
        intron = rc(intron)
    if sum(coverage(intron, units)) / len(intron) < MAX_COVERAGE:
        pool[key] = dict(genome_id=genome_id, seqid=contig, start=start + 1, end=end, strand=strand)

# seeded sample, counted like the telotron loci
ordinary = list(pool.values())
random.Random(SEED).shuffle(ordinary)
sample = ordinary[:SAMPLE_SIZE]
write_tsv(folder / 'ordinary_sample.tsv', sample)
rows = count_junctions(bam, sample, folder / 'ordinary_junctions.tsv')

# sampled introns with >= 1 and >= 3 exact junction reads
summary = dict(
    genome_id=genome_id,
    pool_ordinary_introns=len(ordinary),
    sampled=len(rows),
    ge1=sum(row[4] >= 1 for row in rows),
    ge3=sum(row[4] >= 3 for row in rows),
    seed=SEED,
)
(folder / 'ordinary_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary))
