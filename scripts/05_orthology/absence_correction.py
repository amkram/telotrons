"""Audit every NO_INTRON pair of one genome against comparator annotation and re-call the doubtful
ones (annotation-supported absence correction). Usage: absence_correction.py <stage>.

A pair is reviewed when a comparator CDS model has an annotated intron at the focal position and
phase, or the raw miniprot CIGAR has a target-consuming D, F or G event within 2 focal residues
of the junction (miniprot can encode a short intron as a coding deletion). The qualifying
annotated models decide. If all support continuous CDS, NO_INTRON stays. If all support one
exact-position, same-phase intron and the annotated protein maps back onto the focal intron, the
intron is recovered and classified as TELOTRON, MIXED, ORDINARY or AMBIGUOUS. Anything else
becomes UNRESOLVED_BOUNDARY, which closest.py does not count as eligible.
"""

import hashlib
import json
import os
import re
import shutil
from collections import defaultdict

from Bio.Seq import Seq

from classify_states import TELO_MIN_ARRAY, TELO_PURITY, lineage_motifs
from common import MIN_FLANK_MATCHES, REPO, WORK, dump, header, orthology_dir, parse_stage, read
from common import task_genome, use_scratch, write
from hostlib import open_text
from indexed_reference import Reference
from telotron_ortholog_align import (
    ABSENT_MARGIN,
    FLANK_AA,
    INTRON_TOL,
    MAX_ALT_LOCI,
    parse_miniprot,
    protein_alignment,
    run_miniprot,
)
from call_ortholog_states import GATE_SLOP, MIN_ABSENT_SIM, MIN_SIM, WINDOW_AA, intron_state
from _common import ensure_tool_on_path, find_genome_files, flank_gate, locus_stub, rc, read_fasta

POLICY = 'annotation_supported_absence_correction_v1'
EVENT_WINDOW = 2  # focal residues either side of the junction searched for CIGAR events
TARGET_EVENTS = 'DFG'  # CIGAR operations that consume target DNA without query residues
NEARBY_RESIDUES = 10  # annotated introns mapping this close to the focal residue are recorded
EXCERPT_COLUMNS = 35  # alignment columns either side kept in each decision
DECISION_FIELDS = ['locus', 'ortholog', 'protein', 'old_state', 'new_state', 'reason']


def translate(dna):
    """Translation of the whole codons of a CDS, trailing stops removed."""
    return str(Seq(dna[: len(dna) // 3 * 3]).translate()).rstrip('*')


def paf_records(gff):
    """The ##PAF lines of a miniprot GFF, by query protein."""
    records = defaultdict(list)
    for line in open(gff):
        if line.startswith('##PAF\t'):
            fields = line.rstrip('\n').split('\t')[1:]
            records[fields[0]].append(fields)
    return records


def boundary_events(records, hit, residue, reference):
    """(CIGAR status, indel/frameshift events within EVENT_WINDOW residues of the junction)."""
    matching = [
        p
        for p in records
        if p[5] == hit['contig']
        and p[4] == hit['strand']
        and int(p[2]) + 1 == hit['qmin']
        and int(p[3]) == hit['qmax']
    ]
    if len(matching) != 1:
        return 'ambiguous_PAF', []

    paf = matching[0]
    tags = {tag.split(':', 2)[0]: tag.split(':', 2)[2] for tag in paf[12:]}
    target_start, target_end, strand = int(paf[7]), int(paf[8]), paf[4]
    query_pos = int(paf[2])
    target_used = 0
    events = []
    for length, op in re.findall(r'(\d+)([MIDFGNUV])', tags['cg']):
        length = int(length)
        if op in 'MI':
            query_bp = length
        elif op in 'GUV':
            query_bp = 1
        else:
            query_bp = 0
        if op in 'MD':
            target_bp = 3 * length
        elif op == 'I':
            target_bp = 0
        else:
            target_bp = length

        # an event at the junction, with its genomic interval and DNA
        if op in 'DIFG' and abs(query_pos - residue) <= EVENT_WINDOW:
            if strand == '+':
                start, end = target_start + target_used + 1, target_start + target_used + target_bp
            else:
                start, end = target_end - target_used - target_bp + 1, target_end - target_used
            seq = reference.fetch(paf[5], start - 1, end) if target_bp else ''
            events.append(
                dict(
                    op=op,
                    length=length,
                    target_bp=target_bp,
                    query_boundary=query_pos,
                    delta=query_pos - residue,
                    start=start,
                    end=end,
                    sequence=rc(seq) if strand == '-' else seq,
                )
            )
        query_pos += query_bp
        target_used += target_bp
    assert query_pos == int(paf[3]) and target_used == target_end - target_start, 'CIGAR length'
    return 'checked', events


def annotation_models(gff, surfaces):
    """CDS rows of the comparator annotation on the hit (contig, strand) surfaces, by parent."""
    models = defaultdict(list)
    with open_text(gff) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] != 'CDS' or (fields[0], fields[6]) not in surfaces:
                continue
            attributes = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
            for parent in attributes.get('Parent', attributes.get('ID', '')).split(','):
                models[parent].append((fields, attributes))
    return models


def nearby_annotated_introns(items, seqs, trim, dna, focal, target, hit, residue, phase, reference):
    """Annotated introns of one model that map within NEARBY_RESIDUES of the focal residue."""
    # focal residue count at each target residue's column
    focal_count = 0
    target_count = 0
    target_columns = {}
    for column, (focal_char, target_char) in enumerate(zip(focal, target)):
        if focal_char != '-':
            focal_count += 1
        if target_char != '-':
            target_count += 1
            target_columns[target_count] = (focal_count, column)

    minus = hit['strand'] == '-'
    introns = []
    coding_bp = -trim
    for j in range(len(items) - 1):
        upstream, downstream = items[j][0], items[j + 1][0]
        coding_bp += len(seqs[j])
        target_residue, intron_phase = divmod(coding_bp, 3)
        if minus:
            start, end = int(downstream[4]) + 1, int(upstream[3]) - 1
        else:
            start, end = int(upstream[4]) + 1, int(downstream[3]) - 1
        mapped, column = target_columns.get(target_residue, (None, None))
        if mapped is None or abs(mapped - residue) > NEARBY_RESIDUES:
            continue

        # a gap within 2 columns of the split codon leaves the position uncertain
        near = slice(max(0, column - 2), column + 3)
        gap_near = column is None or '-' in focal[near] or '-' in target[near]
        seq = reference.fetch(hit['contig'], start - 1, end)
        seq = rc(seq) if minus else seq

        # the coding sequence retaining this intron keeps the split-codon frame
        retained = dna[:coding_bp] + seq + dna[coding_bp:]
        introns.append(
            dict(
                start=start,
                end=end,
                bp=end - start + 1,
                sequence=seq,
                splice=seq[:2] + '-' + seq[-2:],
                target_residue=target_residue,
                phase=intron_phase,
                mapped_focal_residue=mapped,
                delta=mapped - residue,
                phase_matches=intron_phase == phase,
                gap_near_boundary=gap_near,
                spliced_stops=translate(dna).count('*'),
                retained_stops=translate(retained).count('*'),
            )
        )
    return introns


def annotated_candidates(models, hit, exon, focal_protein, residue, phase, reference):
    """Annotated CDS models overlapping the hit exon at the junction, translated and aligned."""
    candidates = []
    if not exon:
        return candidates
    for parent, items in models.items():
        overlaps = any(
            fields[0] == hit['contig']
            and fields[6] == hit['strand']
            and int(fields[3]) <= exon[1]
            and int(fields[4]) >= exon[0]
            for fields, _ in items
        )
        if not overlaps:
            continue
        if any(fields[7] not in ['0', '1', '2'] for fields, _ in items):
            candidates.append(
                dict(
                    parent=parent,
                    annotation_status='missing_CDS_phase',
                    phase_continuity_errors=[dict(reason='missing_CDS_phase')],
                    nearby_introns=[],
                )
            )
            continue

        # CDS in transcript order, each phase matching the coding length before it
        items = sorted(items, key=lambda item: int(item[0][3]), reverse=hit['strand'] == '-')
        seqs = [reference.fetch(hit['contig'], int(f[3]) - 1, int(f[4])) for f, _ in items]
        if hit['strand'] == '-':
            seqs = [rc(seq) for seq in seqs]
        coding_bp = -int(items[0][0][7])
        phase_errors = []
        for j, (fields, _) in enumerate(items):
            if j and int(fields[7]) != (-coding_bp) % 3:
                expected = (-coding_bp) % 3
                phase_errors.append(dict(exon=j + 1, observed=int(fields[7]), expected=expected))
            coding_bp += len(seqs[j])

        trim = int(items[0][0][7])
        dna = ''.join(seqs)[trim:]
        protein = translate(dna)
        focal, target = protein_alignment(focal_protein, protein)
        introns = nearby_annotated_introns(
            items, seqs, trim, dna, focal, target, hit, residue, phase, reference
        )
        candidates.append(
            dict(
                parent=parent,
                annotated_protein=items[0][1].get('protein_id', ''),
                phase_continuity_errors=phase_errors,
                nearby_introns=introns,
                cds=[dict(start=int(f[3]), end=int(f[4]), phase=f[7]) for f, _ in items],
                protein_sequence=protein,
                alignment=(focal, target),
            )
        )
    return candidates


def audit_absences(stage, genome_id, absences, map_rows, queries, query_sha):
    """Annotation candidates and CIGAR boundary events of every NO_INTRON pair."""
    source_dir = orthology_dir(stage, genome_id)
    raw = REPO / 'data/raw'
    rows = []
    by_comparator = defaultdict(list)
    for state in absences:
        by_comparator[state['ortholog']].append(state)

    for oid, states in sorted(by_comparator.items()):
        fasta, gff = find_genome_files(oid, str(raw / 'refseq'), str(raw / 'tara'))
        reference = Reference(fasta)
        miniprot_gff = source_dir / '_work' / f'{oid}.{query_sha}.miniprot.gff'
        hits = parse_miniprot(miniprot_gff)
        records = paf_records(miniprot_gff)

        # the hit exon holding each junction
        targets = []
        for state in states:
            map_row = map_rows[(state['telotron'], oid)]
            hit = hits[map_row['protein']]
            residue = int(map_row['intron_residue'])
            exon = next(
                (
                    (start, end, query_start, query_end, identity)
                    for query_start, query_end, start, end, identity in hit['cds_q']
                    if query_start <= residue <= query_end
                ),
                None,
            )
            targets.append((state, map_row, hit, exon))
        surfaces = {(hit['contig'], hit['strand']) for _, _, hit, _ in targets}
        models = annotation_models(gff, surfaces) if gff else {}

        for state, map_row, hit, exon in targets:
            residue, phase = int(map_row['intron_residue']), int(map_row['intron_phase'])
            protein = map_row['protein']
            cigar_status, events = boundary_events(records[protein], hit, residue, reference)
            candidates = annotated_candidates(
                models, hit, exon, queries[protein], residue, phase, reference
            )
            conflicts = [
                intron
                for candidate in candidates
                for intron in candidate['nearby_introns']
                if intron['delta'] == 0
                and intron['phase_matches']
                and not intron['gap_near_boundary']
            ]
            rows.append(
                dict(
                    locus=state['telotron'],
                    ortholog=oid,
                    protein=protein,
                    focal_residue=residue,
                    contig=hit['contig'],
                    strand=hit['strand'],
                    annotation_available=bool(gff),
                    cigar_status=cigar_status,
                    exact_position_phase_conflicts=len(conflicts),
                    boundary_events=events,
                    candidates=candidates,
                )
            )
    return rows


def assess_model(candidate, row, focal_id, reverse_queries):
    """Whether one annotated model supports an intron, continuous CDS, or neither."""
    model = dict(
        parent=candidate['parent'],
        protein=candidate.get('annotated_protein', ''),
        qualified=False,
        reason='invalid_or_missing_CDS_phase',
    )
    if not candidate.get('cds') or candidate.get('phase_continuity_errors'):
        return model

    # the same flank and local-window gates as the state caller
    residue = row['focal_residue']
    focal, target = candidate['alignment']
    columns = [i for i, x in enumerate(focal) if x != '-']
    flanks = [columns[max(0, residue - FLANK_AA) : residue], columns[residue : residue + FLANK_AA]]
    matches = [sum(focal[i] == target[i] and focal[i] not in '-X*' for i in w) for w in flanks]
    ambiguous = sum(focal[i] in 'X*' or target[i] in 'X*' for w in flanks for i in w)
    msa = [(focal_id, focal), ('TARGET', target)]
    window_ok, local_id, n_indel = flank_gate(msa, focal_id, residue, WINDOW_AA, 0, GATE_SLOP)[
        'TARGET'
    ]
    nearby = candidate['nearby_introns']
    model.update(
        matches5=matches[0],
        matches3=matches[1],
        ambiguous_flank_symbols=ambiguous,
        local_gate=[window_ok, local_id, n_indel],
        nearby_introns=nearby,
        reason='fails_local_alignment',
    )
    if min(matches) < MIN_FLANK_MATCHES or ambiguous or not window_ok:
        return model
    model['qualified'] = True

    # an exact-position, same-phase intron, else continuous CDS with no intron nearby
    exact = [
        z for z in nearby if z['delta'] == 0 and z['phase_matches'] and not z['gap_near_boundary']
    ]
    if exact and local_id >= MIN_SIM:
        if len({(z['start'], z['end'], z['phase']) for z in exact}) == 1:
            key = json.dumps([row['ortholog'], candidate['parent'], row['protein']])
            token = 'R' + hashlib.sha256(key.encode()).hexdigest()[:20]
            reverse_queries[token] = candidate['protein_sequence']
            model.update(
                support='intron',
                intron=exact[0],
                reason='exact_position_phase_supported',
                reverse_query=token,
            )
        else:
            model.update(support='uncertain', reason='multiple_exact_boundaries')
    elif not any(abs(z['delta']) <= INTRON_TOL for z in nearby) and local_id >= MIN_ABSENT_SIM:
        # every target-consuming event must lie inside one annotated CDS
        events = [e for e in row['boundary_events'] if e['op'] in TARGET_EVENTS]
        inside = all(
            any(x['start'] <= e['start'] <= e['end'] <= x['end'] for x in candidate['cds'])
            for e in events
        )
        if inside:
            model.update(support='continuous', reason='continuous_annotated_CDS')
        else:
            model.update(support='uncertain', reason='target_event_not_in_one_annotated_CDS')
    else:
        model.update(support='uncertain', reason='nearby_intron_or_low_local_identity')

    center = columns[min(residue, len(columns) - 1)]
    excerpt = slice(max(0, center - EXCERPT_COLUMNS), center + EXCERPT_COLUMNS + 1)
    model.update(alignment_focal=focal[excerpt], alignment_target=target[excerpt])
    return model


def decide(row, models):
    """One reviewed pair: NO_INTRON kept, intron pending reverse mapping, or unresolved."""
    qualified = [m for m in models if m['qualified']]
    supports = {m.get('support') for m in qualified}
    decision = dict(
        locus=row['locus'],
        ortholog=row['ortholog'],
        protein=row['protein'],
        old_state='NO_INTRON',
        new_state='UNRESOLVED_BOUNDARY',
        reason='conflicting_or_uncertain_models' if qualified else 'no_usable_annotation',
        models=models,
        boundary_events=row['boundary_events'],
        annotation_available=row['annotation_available'],
    )
    if supports == {'continuous'}:
        decision.update(new_state='NO_INTRON', reason='independent_continuous_CDS_support')
    elif supports == {'intron'}:
        coordinates = {
            (m['intron']['start'], m['intron']['end'], m['intron']['phase']) for m in qualified
        }
        if len(coordinates) == 1:
            decision.update(
                reason='pending_reverse_mapping', recovered_intron=qualified[0]['intron']
            )
    return decision


def reverse_maps(hit, model, locus):
    """The annotated protein's best hit on the focal genome spans the focal intron cleanly."""
    start, end = int(locus['start']), int(locus['end'])
    residue = model['intron']['target_residue']
    return bool(
        hit
        and hit['contig'] == locus['seqid']
        and hit['strand'] == locus['strand']
        and min(seg[0] for seg in hit['segs']) < start
        and max(seg[1] for seg in hit['segs']) > end
        and hit['n_alternative_loci'] <= MAX_ALT_LOCI
        and not any(
            a <= residue + ABSENT_MARGIN and b >= residue - ABSENT_MARGIN
            for a, b in hit['fs_spans'] + hit['stop_spans']
        )
    )


def finish_recoveries(decisions, reverse_hits, catalogue, motifs):
    """Classify recovered introns that map back to the focal locus, leaving the rest unresolved."""
    recovered = defaultdict(dict)
    for decision in decisions:
        if decision['reason'] != 'pending_reverse_mapping':
            continue
        locus = catalogue[decision['locus']]
        passed = []
        for model in decision['models']:
            if model['qualified']:
                hit = reverse_hits.get(model['reverse_query'])
                model.update(reverse_pass=reverse_maps(hit, model, locus), reverse_hit=hit)
                passed.append(model['reverse_pass'])
        seq = decision['recovered_intron']['sequence']
        if passed and all(passed) and not set(seq) - set('ACGT'):
            state, purity, longest = intron_state(seq, motifs, TELO_PURITY, TELO_MIN_ARRAY)
            decision.update(
                new_state=state,
                reason='annotation_intron_recovered_and_reciprocally_validated',
                repeat_purity=purity,
                longest_array_bp=longest,
                intron_bp=len(seq),
            )
            recovered[decision['locus']][decision['ortholog']] = seq
        elif not all(passed):
            decision['reason'] = 'annotated_protein_reverse_mapping_failed'
        else:
            decision['reason'] = 'ambiguous_bases_in_recovered_intron'
    return recovered


def apply_decisions(states, decisions):
    """Corrected state rows: unresolved and recovered pairs change, all others are copied."""
    by_pair = {(row['telotron'], row['ortholog']): dict(row) for row in states}
    for decision in decisions:
        row = by_pair[(decision['locus'], decision['ortholog'])]
        row['state'] = decision['new_state']
        if decision['new_state'] == 'UNRESOLVED_BOUNDARY':
            row['aligner_call'] = 'UNRESOLVED_BOUNDARY'
        elif decision['new_state'] != 'NO_INTRON':
            local = min(m['local_gate'][1] for m in decision['models'] if m['qualified'])
            row.update(
                aligner_call='ANNOTATION_INTRON',
                repeat_purity=round(decision['repeat_purity'], 4),
                longest_array_bp=decision['longest_array_bp'],
                intron_bp=decision['intron_bp'],
                local_ident=round(local, 4),
            )
    return [by_pair[(row['telotron'], row['ortholog'])] for row in states]


def write_intron_dna(recall_dir, out_dir, recovered):
    """intron_dna_raw/: the recall files, with recovered comparator introns replacing theirs."""
    target = out_dir / 'intron_dna_raw'
    if target.exists():
        shutil.rmtree(target)
    target.mkdir()
    for path in (recall_dir / 'intron_dna_raw').glob('*.fasta'):
        (target / path.name).symlink_to(path)
    for tid, seqs in recovered.items():
        path = target / f'{locus_stub(tid)}.fasta'
        records = list(read_fasta(str(path)).items()) if path.exists() else []
        records = [(name, s) for name, s in records if name.split('|')[0] not in seqs]
        for oid, seq in sorted(seqs.items()):
            records.append((f'{oid}|annotation_recovered_intron|len={len(seq)}', seq))
        path.unlink(missing_ok=True)
        path.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in records))


def main():
    use_scratch()
    ensure_tool_on_path('miniprot')
    stage = parse_stage(__doc__)
    meta = task_genome()
    genome_id = meta['genome_id']
    source_dir = orthology_dir(stage, genome_id)
    recall_dir = orthology_dir(stage, genome_id, 'recall')
    out_dir = orthology_dir(stage, genome_id, 'corrected')
    out_dir.mkdir(parents=True, exist_ok=True)

    # recall states and map rows, mapped catalogue and focal proteins
    states = read(recall_dir / 'states.tsv')
    map_rows = {(r['telotron'], r['ortholog']): r for r in read(recall_dir / 'ortholog_map.tsv')}
    mapped = WORK / stage / 'host' / genome_id / 'mapped_catalogue.tsv'
    catalogue = {row['locus_id']: row for row in read(mapped)}
    query_faa = source_dir / '_work/focal_host_proteins.faa'
    queries = read_fasta(str(query_faa))
    query_sha = hashlib.sha1(query_faa.read_bytes()).hexdigest()[:12]

    # audit every absence and review those with an annotation conflict or target-consuming event
    absences = [row for row in states if row['state'] == 'NO_INTRON']
    audit = audit_absences(stage, genome_id, absences, map_rows, queries, query_sha)
    selected = [
        row
        for row in audit
        if row['exact_position_phase_conflicts']
        or any(event['op'] in TARGET_EVENTS for event in row['boundary_events'])
    ]
    reverse_queries = {}
    decisions = []
    for row in selected:
        models = [assess_model(c, row, genome_id, reverse_queries) for c in row['candidates']]
        decisions.append(decide(row, models))

    # map supporting annotated proteins back onto the focal genome
    reverse_hits = {}
    if reverse_queries:
        validation = out_dir / 'validation'
        validation.mkdir(exist_ok=True)
        query = validation / 'annotated_proteins.faa'
        query.write_text(''.join(f'>{k}\n{v}\n' for k, v in sorted(reverse_queries.items())))
        raw = REPO / 'data/raw'
        fasta, _ = find_genome_files(genome_id, str(raw / 'refseq'), str(raw / 'tara'))
        threads = int(os.environ.get('SLURM_CPUS_PER_TASK', '2'))
        run_miniprot(fasta, str(query), threads, str(validation / 'reverse.miniprot.gff'))
        reverse_hits = parse_miniprot(validation / 'reverse.miniprot.gff')
    motifs = lineage_motifs(genome_id)
    recovered = finish_recoveries(decisions, reverse_hits, catalogue, motifs)

    # corrected states, decisions and intron DNA with the recovered introns
    write(
        out_dir / 'states.tsv',
        apply_decisions(states, decisions),
        header(recall_dir / 'states.tsv'),
    )
    dump(out_dir / 'decisions.json', dict(policy=POLICY, motifs=motifs, decisions=decisions))
    flat = [{field: d[field] for field in DECISION_FIELDS} for d in decisions]
    write(out_dir / 'pair_decisions.tsv', flat, DECISION_FIELDS)
    write_intron_dna(recall_dir, out_dir, recovered)
    changed = sum(d['old_state'] != d['new_state'] for d in decisions)
    print(
        'CORRECTED',
        stage,
        genome_id,
        len(audit),
        'absences',
        len(decisions),
        'reviewed',
        changed,
        'changed',
        flush=True,
    )


if __name__ == '__main__':
    main()
