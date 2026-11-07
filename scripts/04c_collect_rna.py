#!/usr/bin/env python3
"""Join every genome's RNA junction counts to the catalogue, and summarize the ordinary-intron
baseline per genome.

Reads work/04_rna/<genome>/junctions.tsv (04_rna_support.py) and ordinary_summary.json
(04b_ordinary_rna_baseline.py) for each genome of config/rna_runs.tsv. Writes
work/04_rna/RNA_loci.tsv, one row per catalogue locus (status assayed, not_assayed for a locus
without a count in a genome with runs, or no_dataset for a genome without runs), and
work/04_rna/ordinary_baseline.tsv, one row per genome with runs.
"""

import argparse, collections, json, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'lib'))
from audit_its_opportunity import read_tsv, rows_by_genome, write_tsv

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--catalogue', default=str(REPO / 'work/catalogue/catalogue.tsv.gz'))
parser.add_argument('--config', default=str(REPO / 'config/rna_runs.tsv'))
parser.add_argument('--out', default=str(REPO / 'work/04_rna'))
args = parser.parse_args()
out = Path(args.out)
runs = rows_by_genome(args.config)

# exact junction reads per (genome, seqid, start, end)
reads = {}
for genome_id in runs:
    path = out / genome_id / 'junctions.tsv'
    if not path.exists():
        print('no junction counts for', genome_id, '- its loci stay not_assayed', flush=True)
        continue
    for row in read_tsv(path):
        key = (genome_id, row['seqid'], int(row['start']), int(row['end']))
        assert key not in reads, key
        reads[key] = int(row['exact_junction_reads'])

# one row per catalogue locus. A locus without a count gets blank counts, never zero.
loci = []
for locus in read_tsv(args.catalogue):
    genome_id = locus['genome_id']
    count = reads.get((genome_id, locus['seqid'], int(locus['start']), int(locus['end'])))
    if count is not None:
        status = 'assayed'
    elif genome_id in runs:
        status = 'not_assayed'
    else:
        status = 'no_dataset'
    loci.append(
        dict(
            locus_id=locus['locus_id'],
            genome_id=genome_id,
            tier=locus['tier'],
            status=status,
            assayed=int(count is not None),
            exact_junction_reads='' if count is None else count,
            ge1='' if count is None else int(count >= 1),
            ge3='' if count is None else int(count >= 3),
            run_accession=runs[genome_id]['runs'] if genome_id in runs else '',
        )
    )
write_tsv(out / 'RNA_loci.tsv', loci)

# assayed telotrons per genome and tier, and those with >= 1 exact junction read
telotrons = collections.Counter()
for locus in loci:
    if locus['assayed']:
        telotrons[locus['genome_id'], locus['tier'], 'n'] += 1
        telotrons[locus['genome_id'], locus['tier'], 'ge1'] += locus['ge1']

# one baseline row per genome with runs, in config order
baseline = []
for genome_id, row in runs.items():
    ordinary = json.loads((out / genome_id / 'ordinary_summary.json').read_text())
    baseline.append(
        dict(
            genome_id=genome_id,
            organism=row['organism'],
            telotron_primary_ge1=telotrons[genome_id, 'standard', 'ge1'],
            telotron_primary_n=telotrons[genome_id, 'standard', 'n'],
            telotron_additional_ge1=telotrons[genome_id, 'additional', 'ge1'],
            telotron_additional_n=telotrons[genome_id, 'additional', 'n'],
            ordinary_ge1=ordinary['ge1'],
            ordinary_sampled=ordinary['sampled'],
            ordinary_pool=ordinary['pool_ordinary_introns'],
            ordinary_rate=ordinary['ge1'] / ordinary['sampled'],
        )
    )
write_tsv(out / 'ordinary_baseline.tsv', baseline)
print(collections.Counter(locus['status'] for locus in loci))
