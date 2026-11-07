"""Classify comparator introns under the focal lineage's motifs (lib/call_ortholog_states.py)."""

import json
import sys

from common import REPO, WORK, orthology_dir, read, run, sha

CALLER = REPO / 'lib/call_ortholog_states.py'
MOTIF_LENGTHS = range(2, 9)
TELO_PURITY = 0.5  # repeat purity above which an observed intron is TELOTRON
TELO_MIN_ARRAY = 0  # no minimum array length
MIN_FOCAL_PURITY = 0  # no secondary focal purity gate

# motifs without a G or C run of 3 are accepted only for these lineages
SHORT_MOTIF_LINEAGES = {"TTAGG": {"invertebrate"}, "TTAGGC": {"invertebrate"}}


def longest_gc_run(motif):
    """Length of the longest run of G, or of C, in a motif."""
    longest = 0
    for base in "GC":
        run_length = 0
        for x in motif:
            run_length = run_length + 1 if x == base else 0
            longest = max(longest, run_length)
    return longest


def motif_lineage_ok(motif, group):
    """A motif with a G/C run of 3, or a short motif accepted for the genome's group."""
    motif = (motif or "").upper()
    if not motif:
        return False
    if longest_gc_run(motif) >= 3:
        return True
    return group in SHORT_MOTIF_LINEAGES.get(motif, set())


def lineage_motifs(genome_id):
    """The genome's 2-8 bp telomere motifs that are valid for its group (scan_genomes.json)."""
    scan = json.loads((WORK / 'inputs/scan_genomes.json').read_text())
    entry = next(row for row in scan if row['genome_id'] == genome_id)
    motifs = sorted(
        motif
        for motif in entry['motifs']
        if len(motif) in MOTIF_LENGTHS and motif_lineage_ok(motif, entry['group'])
    )
    assert motifs, genome_id
    return motifs


def classify(stage, meta, recall=False):
    """Run the state caller on one genome's orthology (or recall) outputs."""
    genome_id = meta['genome_id']
    directory = orthology_dir(stage, genome_id, 'recall' if recall else '')
    motifs = lineage_motifs(genome_id)

    command = [sys.executable, CALLER, '--ortho-dir', directory, '--focal', genome_id]
    command += ['--motifs', ','.join(motifs), '--min-focal-purity', str(MIN_FOCAL_PURITY)]
    command += ['--telo-purity', str(TELO_PURITY), '--telo-min-array', str(TELO_MIN_ARRAY)]
    command += ['--out', directory / 'states.tsv']
    if recall:
        command.append('--raw-dna')
    run(command, cwd=REPO)
    return dict(
        caller_sha256=sha(CALLER),
        comparator_motifs=motifs,
        comparator_coverage=(
            'maximum coverage by one prespecified focal-lineage motif; no cross-motif union'
        ),
        states=len(read(directory / 'states.tsv')),
        states_sha256=sha(directory / 'states.tsv'),
    )
