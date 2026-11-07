#!/usr/bin/env python3
"""Stage 1: admit annotated introns made of the lineage telomere motif, classify their repeat
architecture, collapse duplicated loci and resolve one strand per locus.
config/strand_corrections.tsv lists curated strand reversals of annotated gene models, applied
after dedup. --collect concatenates the per-genome catalogues into one catalogue.
"""

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'lib'))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    'genomes',
    nargs='*',
    help='genome ids (default: every bearer and control in config/genomes.tsv)',
)
parser.add_argument('--config', default=str(REPO / 'config/genomes.tsv'))
parser.add_argument('--refseq-dir', default=str(REPO / 'data/raw/refseq'))
parser.add_argument('--tara-dir', default=str(REPO / 'data/raw/tara'))
parser.add_argument('--out', default=str(REPO / 'work/01_scan'))
parser.add_argument('--threads', type=int, default=8)
parser.add_argument(
    '--collect',
    action='store_true',
    help='concatenate the per-genome catalogues in --out into --catalogue instead of scanning',
)
parser.add_argument('--catalogue', default=str(REPO / 'work/catalogue/catalogue.tsv.gz'))
args = parser.parse_args()

# keep this run's temporary files under the output directory
OUT = Path(args.out)
SCRATCH = OUT / 'scratch'
OUT.mkdir(parents=True, exist_ok=True)
SCRATCH.mkdir(exist_ok=True)
for variable in ('TMPDIR', 'TMP', 'TEMP'):
    os.environ[variable] = str(SCRATCH)

from _common import find_genome_files, rc, read_tsv, write_tsv
from array_shape import classify
from audit_its_opportunity import fasta
from boundary_policy import apply_boundary_gate
from coverage_policy import apply_coverage_gate
import dedup_telotrons as dedup
import degenerate_linker
from derive_introns import derive
import measure_junctions as junctions
from repeat_policy import coverage, motif_units
from resolve_locus_strands import resolve_strands
from telomere_frac import exact_telomeric_fraction

CORRECTIONS = REPO / 'config/strand_corrections.tsv'
MIN_INTRON_BP = 30
JUNCTION_BP = 12
FLANK_BP = 200


def annotation_support(key, gaps, corrected):
    """Annotated exon pairs bounding an intron, and their basis."""
    support = gaps.get(key, [])
    if support:
        return support, 'same_strand_annotation'

    # a curated strand reversal uses the opposite-strand models, exon lengths swapped
    if corrected:
        seqid, start, end, strand = key
        opposite = gaps.get((seqid, start, end, '-' if strand == '+' else '+'), [])
        if opposite:
            swapped = [
                (model, gene, right, left, kind) for model, gene, left, right, kind in opposite
            ]
            return swapped, 'established_strand_correction'
    return [], 'unsupported'


def find_candidates(introns, contigs, motif, corrections, identity):
    """Introns >= MIN_INTRON_BP over half motif-covered, plus each curated strand reversal."""
    # introns are 0-based half-open, records and locus IDs 1-based closed
    units = tuple(sorted(motif_units(motif)))
    candidates = []
    for seqid, start0, end0, strand in sorted(introns):
        if end0 - start0 < MIN_INTRON_BP:
            continue
        intron = contigs[seqid][start0:end0]
        if strand == '-':
            intron = rc(intron)
        if sum(coverage(intron, units)) / len(intron) <= 0.5:
            continue

        # a corrected locus also gets a record on its curated strand
        locus_id = f"{identity['genome_id']}|{seqid}:{start0 + 1}-{end0}"
        strands = [strand]
        curated = corrections.get(locus_id)
        if curated and curated != strand and (seqid, start0, end0, curated) not in introns:
            strands.append(curated)
        for record_strand in strands:
            candidates.append(
                dict(
                    identity,
                    locus_id=locus_id,
                    seqid=seqid,
                    start=start0 + 1,
                    end=end0,
                    strand=record_strand,
                )
            )
    return candidates


def measure_locus(candidate, motif, contigs, gaps, corrected):
    """Shape, junction copy core, gates and admission filters of one candidate, or None."""
    key = (candidate['seqid'], candidate['start'], candidate['end'], candidate['strand'])
    flanked = junctions.window(contigs, key, FLANK_BP)
    if flanked is None:
        return None
    left, intron, right = flanked

    # junction copy core on flanks trimmed to the supporting exons
    support, basis = annotation_support(key, gaps, corrected)
    exon_flanked = junctions.bound_flanks(junctions.window(contigs, key), support)
    core = junctions.measure(*exon_flanked, motif)

    # array shape, then the boundary and coverage gates on the architecture
    shape = classify(intron, motif)
    gate = apply_boundary_gate(
        shape['architecture'], len(intron), shape['arrays'], core, bool(support)
    )
    coverage_gate = apply_coverage_gate(gate['architecture'], shape['covered_bp'], len(intron))

    row = dict(candidate)
    row.update(
        motif=motif,
        tier=shape['tier'],
        covered_bp=shape['covered_bp'],
        repeat_fraction=shape['repeat_fraction'],
        shape_architecture=shape['architecture'],
        shape_reason=shape['reason'],
        arrays_json=json.dumps(shape['arrays'], separators=(',', ':')),
        n_arrays=len(shape['arrays']),
        unexplained_repeat_bp=shape['unexplained_repeat_bp'],
        repeat_bp_outside_arrays=shape['repeat_bp_outside_arrays'],
        linker_start=shape['linker_start'],
        linker_end=shape['linker_end'],
        longest_nonrepeat_linker=shape['longest_nonrepeat_linker'],
        ignored_edge_repeat_bp=shape['ignored_edge_repeat_bp'],
        bridged_gap_bp=shape['bridged_gap_bp'],
    )
    row.update(gate)
    row.update({'tsd_' + name: value for name, value in core.items()})
    row['tsd_annotation_support_basis'] = basis
    row['tsd_boundary_censored'] = int(bool(core['core_censored'] and gate['tsd_adjustment_used']))
    row.update(coverage_gate)

    # degenerate-linker screen, run only on linker architectures
    reasons = [shape['reason'], gate['boundary_gate_reason'], coverage_gate['coverage_gate_reason']]
    row.update(
        degenerate_linker_screen='not_applicable',
        degenerate_residual_start0='',
        degenerate_residual_end0='',
        degenerate_residual_core_bp='',
    )
    if '-L-' in row['architecture']:
        screen = degenerate_linker.evaluate(intron, motif, *shape['arrays'])
        row.update(
            degenerate_linker_screen='retained' if screen['linker_retained'] else 'excluded',
            degenerate_residual_start0=screen['left_extended_to0'],
            degenerate_residual_end0=screen['right_extended_to0'],
            degenerate_residual_core_bp=screen['residual_nonrepeat_bp'],
        )
        if not screen['linker_retained']:
            row['architecture'] = 'mixed_unclear'
            reasons.append('degenerate_repeat_extension_leaves_lt18bp_core')
    row['architecture_reason'] = ';'.join(reason for reason in reasons if reason)

    # junction (12 bp) and flank (200 bp) repeat content
    units = motif_units(motif)
    flank_repeat_bp = sum(coverage(left, units)) + sum(coverage(right, units))
    row.update(
        flank_cont5=exact_telomeric_fraction(left[-JUNCTION_BP:], motif),
        flank_cont3=exact_telomeric_fraction(right[:JUNCTION_BP], motif),
        flank_repeat_fraction=flank_repeat_bp / (len(left) + len(right)) if left or right else 1.0,
    )

    # admission filters
    reject = []
    if min(len(left), len(right)) < JUNCTION_BP:
        reject.append('insufficient_junction_flank')
    if max(row['flank_cont5'], row['flank_cont3']) >= 0.5:
        reject.append('repeat_continues_across_junction')
    if row['flank_repeat_fraction'] > 0.2:
        reject.append('flanking_repeat_fraction_gt20')
    if row['repeat_fraction'] <= 0.5:
        reject.append('coverage_le50')
    row.update(
        filter_pass=int(not reject),
        filter_reasons=';'.join(reject),
        splice_class=intron[:2] + '-' + intron[-2:],
        intron_len=len(intron),
        intron_sha256=hashlib.sha256(intron.encode()).hexdigest(),
    )
    return row


def scan_genome(genome_id, meta, corrections):
    """Admit, deduplicate and strand-resolve one genome's loci, then write its tables and stats."""
    started = time.time()
    motif = meta['motif']
    fasta_path, gff = find_genome_files(genome_id, args.refseq_dir, args.tara_dir)
    assert fasta_path and gff, (genome_id, fasta_path, gff)

    # candidate introns from the annotation
    introns = derive(gff)
    wanted = {seqid for seqid, _, _, _ in introns}
    contigs = {seqid: seq.upper() for seqid, seq in fasta(fasta_path) if seqid in wanted}
    gaps = junctions.annotation_gaps(gff)
    identity = dict(genome_id=genome_id, organism=meta['label'], clade=meta['lineage'])
    candidates = find_candidates(introns, contigs, motif, corrections, identity)

    # measure each candidate and keep those passing every admission filter
    admitted = []
    dropped = Counter()
    for candidate in candidates:
        corrected = corrections.get(candidate['locus_id']) == candidate['strand']
        row = measure_locus(candidate, motif, contigs, gaps, corrected)
        if row is None:
            dropped['window_unavailable'] += 1
        elif row['filter_pass']:
            admitted.append(row)
        else:
            # rejected loci are tallied under their first reason
            dropped[row['filter_reasons'].split(';')[0]] += 1

    # group duplicated loci and keep one representative per group
    groups = dedup.dedup_groups(genome_id, admitted, contigs, SCRATCH, args.threads)
    representatives = dedup.choose_representatives(admitted, groups)
    for i, row in enumerate(admitted):
        row.update(
            dedup_keep=int(i == representatives[i]),
            dedup_group=groups[i],
            representative_locus=admitted[representatives[i]]['locus_id'],
        )

    # one record per locus, with strand '.' where its strand records conflict
    retained = [row for row in admitted if row['dedup_keep']]
    resolved = resolve_strands(retained, corrections)

    stats = dict(
        genome_id=genome_id,
        motif=motif,
        candidates=len(candidates),
        filter_pass=len(admitted),
        dropped=dict(dropped),
        after_dedup=len(retained),
        after_strand=len(resolved),
        tiers=dict(Counter(row['tier'] for row in resolved)),
        seconds=round(time.time() - started, 1),
    )
    write_tsv(OUT / f'{genome_id}.admitted.tsv.gz', admitted)
    write_tsv(OUT / f'{genome_id}.catalogue.tsv.gz', resolved)
    (OUT / f'{genome_id}.stats.json').write_text(json.dumps(stats, indent=1))
    print('DONE', json.dumps(stats), flush=True)


def collect(genome_ids):
    """Concatenate the per-genome catalogues. Every genome must have been scanned."""
    rows = []
    for genome_id in genome_ids:
        stats = json.loads((OUT / f'{genome_id}.stats.json').read_text())
        path = OUT / f'{genome_id}.catalogue.tsv.gz'
        genome_rows = read_tsv(path) if path.exists() else []
        assert len(genome_rows) == stats['after_strand'], genome_id
        rows.extend(genome_rows)
    Path(args.catalogue).parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.catalogue, rows)
    print(f'{len(rows)} loci from {len(genome_ids)} genomes -> {args.catalogue}', flush=True)


if __name__ == '__main__':
    genomes = {row['genome_id']: row for row in read_tsv(args.config)}
    genome_ids = args.genomes or [
        genome_id for genome_id, meta in genomes.items() if meta['role'] in ('bearer', 'control')
    ]
    if args.collect:
        collect(genome_ids)
    else:
        corrections = {row['locus_id']: row['strand'] for row in read_tsv(CORRECTIONS)}
        for genome_id in genome_ids:
            scan_genome(genome_id, genomes[genome_id], corrections)
