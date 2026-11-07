"""Screen each control's annotated introns (>=30 bp) for lineage-motif coverage and apply the
junction and flank filters to those >50% covered.
"""

import argparse
import json
from pathlib import Path
import sys
import time

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'lib'))
from _common import find_genome_files, read_tsv
from audit_its_opportunity import fasta
import measure_junctions as junctions
from repeat_policy import coverage, motif_units
from telomere_frac import exact_telomeric_fraction

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('genomes', nargs='+')
parser.add_argument('--config', default=str(REPO / 'config/genomes.tsv'))
parser.add_argument('--refseq-dir', default=str(REPO / 'data/raw/refseq'))
parser.add_argument('--tara-dir', default=str(REPO / 'data/raw/tara'))
parser.add_argument('--out', default=str(REPO / 'work/02_control_screen'))
args = parser.parse_args()

OUT = Path(args.out)
OUT.mkdir(parents=True, exist_ok=True)
genomes = {row['genome_id']: row for row in read_tsv(args.config)}


def junction_filters(left, right, motif, units):
    """Names of the junction and flank filters an intron fails."""
    junction5 = exact_telomeric_fraction(left[-12:], motif)
    junction3 = exact_telomeric_fraction(right[:12], motif)
    flank_repeat_bp = sum(coverage(left, units)) + sum(coverage(right, units))
    flank_fraction = flank_repeat_bp / max(1, len(left) + len(right))
    checks = (
        ('insufficient_junction_flank', min(len(left), len(right)) < 12),
        ('repeat_continues_across_junction', max(junction5, junction3) >= 0.5),
        ('flanking_repeat_fraction_gt20', flank_fraction > 0.2),
    )
    return [name for name, failed in checks if failed]


def scan(genome_id):
    """Coverage of every annotated intron of one control, and the filters of those >50% covered."""
    started = time.time()
    motif = genomes[genome_id]['motif']
    units = tuple(sorted(motif_units(motif)))

    # report missing annotation or assembly files instead of failing
    fasta_path, gff = find_genome_files(genome_id, args.refseq_dir, args.tara_dir)
    if not fasta_path or not gff:
        return dict(
            genome_id=genome_id,
            status='files_missing',
            fa=int(bool(fasta_path)),
            gff=int(bool(gff)),
        )
    gaps = junctions.annotation_gaps(gff)
    contigs = {seqid: seq.upper() for seqid, seq in fasta(fasta_path)}

    # motif coverage of every intron of at least 30 bp (gap keys are 1-based closed)
    scanned = 0
    short = 0
    max_fraction = 0
    hits = []
    for seqid, start, end, strand in gaps:
        if end - start + 1 < 30:
            short += 1
            continue
        flanked = junctions.window(contigs, (seqid, start, end, strand), 200)
        if flanked is None:
            continue
        left, intron, right = flanked
        scanned += 1
        fraction = sum(coverage(intron, units)) / len(intron)
        max_fraction = max(max_fraction, fraction)

        # introns over 50% covered get the junction and flank filters
        if fraction > 0.5:
            reasons = junction_filters(left, right, motif, units)
            hits.append(
                dict(
                    locus=f'{genome_id}|{seqid}:{start}-{end}',
                    len=len(intron),
                    fraction=round(fraction, 3),
                    pass_filters=int(not reasons),
                    reasons=';'.join(reasons),
                )
            )

    return dict(
        genome_id=genome_id,
        status='scanned',
        motif=motif,
        annotated_unique_introns=scanned,
        shorter_than_30=short,
        coverage_gt50=len(hits),
        pass_filters=sum(hit['pass_filters'] for hit in hits),
        max_fraction=round(max_fraction, 3),
        hits=hits,
        seconds=round(time.time() - started, 1),
    )


if __name__ == '__main__':
    # print a summary per genome and keep the per-intron hits in the JSON file
    for genome_id in args.genomes:
        result = scan(genome_id)
        print(
            json.dumps({key: value for key, value in result.items() if key != 'hits'}), flush=True
        )
        (OUT / f'{genome_id}.rescan.json').write_text(json.dumps(result, indent=1))
