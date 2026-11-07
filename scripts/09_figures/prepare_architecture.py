"""Architecture counts by lineage and tier."""

from pathlib import Path
import collections
import csv
import gzip
import hashlib
import json
import os

REPO = Path(__file__).resolve().parents[2]
PACKAGE = Path(os.environ['TELO_PACKAGE'])

ARCH = ['R', 'F', 'R-F', 'F-R', 'R-L-F', 'F-L-R']
CLADES = [
    'haptophyte',
    'cryptophyte',
    'dinoflagellate',
    'chrysophyte',
    'apicomplexa',
    'nematode',
    'stramenopile',
]


def read(path):
    with gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(name, rows):
    with open(PACKAGE / 'tables' / name, 'w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


def architecture_counts(loci):
    tally = collections.Counter(locus['architecture'] for locus in loci)
    return dict(
        loci=len(loci),
        genomes=len({locus['genome_id'] for locus in loci}),
        clean=sum(tally[architecture] for architecture in ARCH),
        **{architecture: tally[architecture] for architecture in ARCH},
        mixed_unclear=tally['mixed_unclear'],
    )


catalogue = read(REPO / 'work/catalogue/catalogue.tsv.gz')
by_id = {locus['locus_id']: locus for locus in catalogue}
primary = [locus for locus in catalogue if locus['tier'] == 'standard']
additional = [locus for locus in catalogue if locus['tier'] == 'additional']

# every classified locus must meet the architecture rules: outer arrays within 8 bp of the
# boundaries, >85% coverage (>50% with a linker) and a linker core of at least 18 bp
classified = [locus for locus in catalogue if locus['architecture'] in ARCH]
for locus in classified:
    linker = '-L-' in locus['architecture']
    assert float(locus['donor_boundary_distance_bp']) <= 8
    assert float(locus['acceptor_boundary_distance_bp']) <= 8
    assert float(locus['repeat_fraction']) > (0.5 if linker else 0.85)
    if linker:
        assert int(locus['degenerate_residual_core_bp']) >= 18

# Figure 2A: count of each architecture among each lineage's classified loci
rows = []
for architecture in ARCH:
    for clade in CLADES:
        lineage = [locus for locus in classified if locus['clade'] == clade]
        count = sum(locus['architecture'] == architecture for locus in lineage)
        rows.append(
            dict(
                clade=clade,
                architecture=architecture,
                count=count,
                clean_denominator=len(lineage),
                percent=100 * count / len(lineage) if lineage else 0,
            )
        )
write('architecture_counts.tsv', rows)

# architecture counts by coverage tier, overall and per lineage
tiers = []
clades = []
for tier, loci in [('gt85', primary), ('gt50_to85', additional), ('gt50', catalogue)]:
    tiers.append(dict(tier=tier, **architecture_counts(loci)))
    for clade in CLADES:
        lineage = [locus for locus in loci if locus['clade'] == clade]
        clades.append(dict(clade=clade, tier=tier, **architecture_counts(lineage)))
write('tier_summary.tsv', tiers)
write('clade_summary.tsv', clades)

# the curated Figure 2B examples must match their catalogue loci
examples = json.loads((PACKAGE / 'inputs/examples_matched_verified.json').read_text())
for example in examples:
    locus = by_id[example['pid']]
    assert locus['architecture'] == example['architecture']
    assert sha256(example['context']['intron']) == locus['intron_sha256']
    for key in ['motif', 'strand', 'repeat_fraction', 'tsd_core_score', 'tsd_core_sequence']:
        assert locus[key] == example['catalogue'][key], (example['pid'], key)

print(
    json.dumps(
        dict(
            primary=len(primary),
            additional=len(additional),
            classified_primary=architecture_counts(primary)['clean'],
            classified_additional=architecture_counts(additional)['clean'],
        )
    )
)
