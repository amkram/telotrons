#!/usr/bin/env python3
"""ITS density of every genome of config/genomes.tsv, and of each bearer-control pair.

Density is ITS components per Mb of eligible A/C/G/T sequence outside annotated genes, terminal
telomeric caps and projected gene spans, from projection.py's work/07_its/proj/<genome>/summary.tsv.
Writes work/07_its/figure_points.tsv and pairs.tsv.
"""

import argparse, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'lib'))
from audit_its_opportunity import read_tsv, rows_by_genome, write_tsv

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--config', default=str(REPO / 'config/genomes.tsv'))
parser.add_argument('--out', default=str(REPO / 'work/07_its'))
args = parser.parse_args()
out = Path(args.out)
genomes = rows_by_genome(args.config)

# one point per genome, in config order
points = []
for genome_id, genome in genomes.items():
    [row] = read_tsv(out / 'proj' / genome_id / 'summary.tsv')
    n_components = int(row['n_components'])
    eligible = int(row['eligible_acgt_bp'])
    points.append(
        dict(
            genome_id=genome_id,
            role=genome['role'],
            label=genome['label'],
            clade=genome['lineage'],
            assembly_class=genome['assembly_class'],
            ITS_components=n_components,
            eligible_ACGT_bp=eligible,
            density=n_components * 1e6 / eligible,
            pair_bearer_id=genome['pair_bearer_id'],
            pair_control_id=genome['pair_control_id'],
        )
    )
write_tsv(out / 'figure_points.tsv', points)

# one pair per bearer, grouped by the bearer's lineage (each control serves one bearer)
density = {point['genome_id']: point['density'] for point in points}
bearers = [genome_id for genome_id, genome in genomes.items() if genome['role'] == 'bearer']
controls = [genome_id for genome_id, genome in genomes.items() if genome['role'] == 'control']
assert sorted(genomes[bearer]['pair_control_id'] for bearer in bearers) == sorted(controls)
pairs = []
for bearer in bearers:
    control = genomes[bearer]['pair_control_id']
    assert genomes[bearer]['pair_bearer_id'] == genomes[control]['pair_bearer_id'] == bearer
    pairs.append(
        dict(
            group=genomes[bearer]['lineage'],
            bearer_id=bearer,
            control_id=control,
            bearer_density=density[bearer],
            control_density=density[control],
            density_difference=density[bearer] - density[control],
        )
    )
write_tsv(out / 'pairs.tsv', pairs)
smallest = min(pair['density_difference'] for pair in pairs)
print(len(points), 'genomes;', len(pairs), 'pairs; smallest difference', smallest)
