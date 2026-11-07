"""Audit intron position and phase against the annotated models, under two gap penalties."""

import json
import os
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

from exact_position import OUT, read_json_gz, sha, write_json_gz

warnings.filterwarnings('ignore', message='Partial codon')
WORKERS = int(os.environ.get('AUDIT_WORKERS', '4'))
FLANK_AA = 20  # focal residues either side whose matches are counted
MIN_FLANK_MATCHES = 12  # per flank
GATE_AA = 10  # residues either side of the junction without disallowed indels
BOUNDARY_SLOP = 2  # alignment columns either side of the junction where an indel is allowed
EXCERPT_COLUMNS = 45  # alignment columns either side kept for review
NEAR_SITE_RESIDUES = 12  # model introns mapping this close to the focal codon are reported

# two global BLOSUM62 aligners differing only in gap-open penalty
ALIGNERS = []
for gap_open in (-11, -9):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = -1
    aligner.mode = 'global'
    ALIGNERS.append(aligner)


@lru_cache(maxsize=512)
def align(focal, target):
    """Top (focal, target) alignment strings under each gap penalty."""
    if focal == target:
        return ((focal, target), (focal, target))
    alignments = [aligner.align(focal, target)[0] for aligner in ALIGNERS]
    return tuple((str(alignment[0]), str(alignment[1])) for alignment in alignments)


def translate(dna):
    """Translation of the whole codons of a CDS, trailing stops removed."""
    return str(Seq(dna[: len(dna) // 3 * 3]).translate()).rstrip('*')


def local_check(focal, target, residue):
    """Flank matches, ambiguity and indels around the junction in one alignment."""
    columns = [i for i, x in enumerate(focal) if x != '-']
    flanks = [columns[max(0, residue - FLANK_AA) : residue], columns[residue : residue + FLANK_AA]]
    matches = [sum(focal[i] == target[i] and focal[i] not in '-X*' for i in w) for w in flanks]
    ambiguous = sum(focal[i] in 'X*' or target[i] in 'X*' for w in flanks for i in w)
    if not 0 < residue < len(columns):
        return dict(qualified=False, matches=matches, ambiguous=ambiguous)

    # indels within GATE_AA residues are disallowed, except within BOUNDARY_SLOP columns
    # of the junction
    center = columns[residue]
    low = columns[max(0, residue - GATE_AA)]
    high = columns[min(len(columns) - 1, residue + GATE_AA - 1)]
    gaps = [
        i
        for i in range(low, high + 1)
        if (focal[i] == '-' or target[i] == '-') and abs(i - center) > BOUNDARY_SLOP
    ]

    excerpt_start = max(0, center - EXCERPT_COLUMNS)
    return dict(
        qualified=min(matches) >= MIN_FLANK_MATCHES and ambiguous == 0 and not gaps,
        matches=matches,
        ambiguous=ambiguous,
        disallowed_gap_columns=[i - center for i in gaps],
        focal=focal[excerpt_start : center + EXCERPT_COLUMNS + 1],
        target=target[excerpt_start : center + EXCERPT_COLUMNS + 1],
        boundary_column=min(EXCERPT_COLUMNS, center),
        focal_start_residue=sum(x != '-' for x in focal[:excerpt_start]),
        target_start_residue=sum(x != '-' for x in target[:excerpt_start]),
    )


def mapped_site(focal, target, target_residue, target_phase, residue, phase):
    """Where a target intron site (codon, phase) lands on the focal protein in one alignment."""
    # alignment column of each target residue, and focal residue index per column
    target_columns = [i for i, x in enumerate(target) if x != '-']
    focal_index = {}
    n = 0
    for i, x in enumerate(focal):
        if x != '-':
            focal_index[i] = n
            n += 1
    if not 0 <= target_residue < len(target_columns):
        return dict(exact=False, strict=False, reason='outside_alignment')

    column = target_columns[target_residue]
    mapped = focal_index.get(column)
    direct = mapped is not None

    # a phase-0 site sits between residues: the previous one must align directly too
    if target_phase == 0:
        previous = target_columns[target_residue - 1] if target_residue > 0 else None
        direct = (
            direct
            and target_residue > 0
            and previous + 1 == column
            and focal_index.get(previous) == mapped - 1
        )
    near = slice(max(0, column - 2), column + 3)
    gap = '-' in focal[near] or '-' in target[near]
    exact = direct and mapped == residue and target_phase == phase
    return dict(
        mapped_focal_codon=mapped,
        delta=None if mapped is None else mapped - residue,
        target_phase=target_phase,
        phase_matches=target_phase == phase,
        direct_site_mappable=bool(direct),
        gap_within2=gap,
        exact=bool(exact),
        strict=bool(exact and not gap),
    )


def reading_frame(cds, focal_protein):
    """(bases trimmed before the first codon, fields to record), or None when unresolvable."""
    raw = ''.join(c['sequence'] for c in cds)
    if all(c['phase'] in ('0', '1', '2') for c in cds):
        return int(cds[0]['phase']), dict(phase_source='annotation')

    # missing phase: the unique best-scoring frame against the focal protein
    candidates = []
    for frame in range(3):
        protein = translate(raw[frame:])
        if protein:
            candidates.append((float(ALIGNERS[0].score(focal_protein, protein)), frame))
    candidates.sort(reverse=True)
    if not candidates or len(candidates) > 1 and candidates[0][0] == candidates[1][0]:
        return None
    fields = dict(phase_source='inferred_unique_best_protein_frame', frame_scores=candidates)
    return candidates[0][1], fields


def check_model(request, model):
    """Translation, phase consistency, local alignment checks and mapped introns of one model."""
    cds = model['cds']
    record = {k: model[k] for k in ['parent', 'protein_id', 'product']}
    frame = reading_frame(cds, request['focal_protein'])
    if frame is None:
        return dict(record, status='frame_unresolved', qualified=False)
    trim, fields = frame
    record.update(fields)

    # each exon's annotated phase must match the CDS length before it
    coding_bp = -trim
    errors = []
    for i, c in enumerate(cds):
        if i and c['phase'] in ('0', '1', '2') and int(c['phase']) != (-coding_bp) % 3:
            errors.append(dict(exon=i, observed=int(c['phase']), expected=(-coding_bp) % 3))
        coding_bp += len(c['sequence'])

    # qualified: both alignments pass the local check and no phase errors
    protein = translate(''.join(c['sequence'] for c in cds)[trim:])
    record.update(
        trim=trim,
        phase_errors=errors,
        protein_length=len(protein),
        internal_stops=protein.count('*'),
    )
    if not protein:
        return dict(record, status='empty_protein', qualified=False)
    residue, phase = request['focal_residue'], request['focal_phase']
    alignments = align(request['focal_protein'], protein)
    checks = [local_check(focal, target, residue) for focal, target in alignments]
    record.update(
        qualified=all(check['qualified'] for check in checks) and not errors,
        local_checks=checks,
        status='evaluated',
    )

    # map each intron to a target codon and phase, keeping the stored intron and those
    # within NEAR_SITE_RESIDUES of the focal codon
    introns = []
    coding_bp = -trim
    for i, intron in enumerate(model['introns']):
        coding_bp += len(cds[i]['sequence'])
        target_residue, target_phase = divmod(coding_bp, 3)
        sites = [
            mapped_site(focal, target, target_residue, target_phase, residue, phase)
            for focal, target in alignments
        ]
        stored = intron['start'] == request['start'] and intron['end'] == request['end']
        near = any(
            site.get('delta') is not None and abs(site['delta']) <= NEAR_SITE_RESIDUES
            for site in sites
        )
        if stored or near:
            introns.append(
                dict(
                    start=intron['start'],
                    end=intron['end'],
                    bp=intron['bp'],
                    stored_boundary=stored,
                    target_codon=target_residue,
                    target_phase=target_phase,
                    sites=sites,
                    exact_both=all(site['exact'] for site in sites),
                    strict_both=all(site['strict'] for site in sites),
                )
            )
    record['introns'] = introns
    return record


def site_status(request, models):
    """Status of one site from its qualified models (first matching rule wins)."""
    qualified = [m for m in models if m['qualified']]
    stored = [(m, z) for m in qualified for z in m['introns'] if z['stored_boundary']]
    strict = [(m, z) for m, z in stored if z['strict_both']]
    direct = [(m, z) for m, z in stored if z['exact_both']]
    other_exact = [
        (m, z)
        for m in qualified
        for z in m['introns']
        if z['strict_both'] and not z['stored_boundary']
    ]
    disagreement = [(m, z) for m, z in stored if not z['exact_both']]

    if strict and disagreement:
        return 'conflicting_qualified_models', strict, qualified
    if strict:
        return 'exact_position_phase', strict, qualified
    if direct:
        return 'exact_direct_site_near_alignment_gap', strict, qualified
    if other_exact:
        return 'exact_site_different_genomic_boundary', strict, qualified
    if stored:
        # stable: every site mapped directly, away from gaps, to one (delta, phase)
        stable = [
            z
            for m, z in stored
            if all(v.get('direct_site_mappable') and not v.get('gap_within2') for v in z['sites'])
            and len({(v.get('delta'), v.get('target_phase')) for v in z['sites']}) == 1
        ]
        if stable:
            return 'different_position_or_phase', strict, qualified
        return 'position_alignment_uncertain', strict, qualified
    if not request['annotation_available']:
        return 'annotation_unavailable', strict, qualified
    if not models:
        return 'no_overlapping_coding_model', strict, qualified
    if not qualified:
        return 'annotation_alignment_or_phase_unqualified', strict, qualified
    return 'stored_intron_not_in_qualified_annotation', strict, qualified


def analyze(request):
    """Check every overlapping model of one site, then summarize them into one status."""
    models = [check_model(request, model) for model in request['models']]
    status, strict, qualified = site_status(request, models)
    expected = request['expected_intron']
    return dict(
        key=request['key'],
        locus_id=request['locus_id'],
        genome_id=request['genome_id'],
        role=request['role'],
        status=status,
        original_annotated_phase=any(m['phase_source'] == 'annotation' for m, z in strict),
        inferred_phase_only=(
            bool(strict) and all(m['phase_source'] != 'annotation' for m, z in strict)
        ),
        genomic_sequence_matches=(
            request['genomic_sequence'].upper() == expected.upper() if expected else None
        ),
        qualified_models=len(qualified),
        models=models,
    )


def attach_junctions(requests, loci):
    """Add the focal junction (residue, phase) of its locus to every request."""
    by_locus = {locus['locus_id']: locus for locus in loci}
    for request in requests:
        locus = by_locus[request['locus_id']]
        request.update(
            focal_residue=int(locus['intron_residue']), focal_phase=int(locus['intron_phase'])
        )


def main():
    # inputs must match the annotation extraction certificate
    evidence = read_json_gz('evidence.json.gz')
    requests = read_json_gz('annotations.json.gz')['rows']
    certificate = json.loads((OUT / 'annotation_extraction_complete.json').read_text())
    assert sha(OUT / 'annotations.json.gz') == certificate['annotations_sha256']
    assert sha(OUT / 'evidence.json.gz') == certificate['evidence_sha256']
    attach_junctions(requests, evidence['loci'])

    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        for result in pool.map(analyze, requests, chunksize=4):
            results.append(result)
            if len(results) % 50 == 0:
                print('MODELS', len(results), '/', len(requests), flush=True)

    # every stored intron must match the genome
    assert all(result['genomic_sequence_matches'] for result in results)
    write_json_gz('model_review.json.gz', dict(rows=results))

    inputs = (OUT / 'evidence.json.gz', OUT / 'annotations.json.gz')
    complete = dict(
        complete=True,
        requests=len(results),
        statuses=dict(Counter(result['status'] for result in results)),
        input_sha256={str(path): sha(path) for path in inputs},
        output_sha256=sha(OUT / 'model_review.json.gz'),
        script_sha256=sha(__file__),
    )
    (OUT / 'model_review_complete.json').write_text(json.dumps(complete, indent=2) + '\n')
    print(json.dumps(complete, indent=2), flush=True)


if __name__ == '__main__':
    main()
