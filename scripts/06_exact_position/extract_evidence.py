"""Collect the evidence of every closest-ORDINARY primary locus (read-only: calls never change)."""

import hashlib
import json
from collections import Counter, defaultdict

from exact_position import (
    ELIGIBLE_STATES,
    ORTHOLOGY,
    OUT,
    RESOLVED_STATES,
    RNA_LOCI,
    check_scratch,
    primary_closest,
    read_tsv,
    sha,
    write_json_gz,
)
from _common import locus_stub, read_aln
from telotron_ortholog_align import INTRON_TOL, parse_miniprot

ALIGNMENT_WINDOW_AA = 35  # focal residues either side of the junction kept for review
FLANK_AA = 20  # residues either side over which flank matches are recounted


def read_alignment(path):
    """Genome ID (first header field) to sequence of an alignment or FASTA. Empty if missing."""
    return {name.split('|')[0]: seq for name, seq in read_aln(path)}


def tied_ordinary_comparators(states):
    """Comparators tied at the best eligible identity, all of which must be ORDINARY."""
    eligible = [state for state in states if state['state'] in ELIGIBLE_STATES]
    best = max(float(state['protein_ident']) for state in eligible)
    selected = {state['ortholog'] for state in eligible if float(state['protein_ident']) == best}
    assert all(state['state'] == 'ORDINARY' for state in eligible if state['ortholog'] in selected)
    return sorted(selected)


def nearby_prediction(hit, residue):
    """The hit's predicted introns within INTRON_TOL residues of the junction, and the longest."""
    nearby = [
        dict(
            query_residue=position,
            delta=position - residue,
            start=start,
            end=end,
            bp=end - start + 1,
        )
        for position, start, end in hit['introns']
        if abs(position - residue) <= INTRON_TOL
    ]
    return nearby, max(nearby, key=lambda intron: intron['bp']) if nearby else None


def locus_pairs(tid, genome_id, states, map_rows, corrections, parsed_gffs, source_hashes):
    """(protein, residue, intron DNA, pairs), one pair per tied comparator with its evidence."""
    recall_dir = ORTHOLOGY / 'recall/primary/orthology' / genome_id
    corrected_dir = ORTHOLOGY / 'corrected/primary/orthology' / genome_id
    raw_path = corrected_dir / 'intron_dna_raw' / f'{locus_stub(tid)}.fasta'
    intron_seqs = read_alignment(raw_path)
    if raw_path.exists():
        source_hashes[str(raw_path)] = sha(raw_path)

    # the focal protein alignment, with a window of ALIGNMENT_WINDOW_AA residues each side
    protein = map_rows[0]['protein'] if map_rows else ''
    msa_path = recall_dir / 'protein_msa' / f'{genome_id}__{protein}.aln.fasta'
    msa = read_alignment(msa_path)
    if msa_path.exists():
        source_hashes[str(msa_path)] = sha(msa_path)
    focal = msa.get(genome_id, '')
    residue = int(map_rows[0]['intron_residue']) if map_rows else 0
    columns = [i for i, residue_char in enumerate(focal) if residue_char != '-']
    center = columns[residue] if residue < len(columns) else len(focal)
    low = columns[max(0, residue - ALIGNMENT_WINDOW_AA)] if columns else 0
    high = columns[min(len(columns) - 1, residue + ALIGNMENT_WINDOW_AA)] + 1 if columns else 0
    flanks = [columns[max(0, residue - FLANK_AA) : residue], columns[residue : residue + FLANK_AA]]

    pairs = []
    for oid in tied_ordinary_comparators(states):
        map_row = next(row for row in map_rows if row['ortholog'] == oid)
        state = next((row for row in states if row['ortholog'] == oid), None)
        target = msa.get(oid, '')

        matches = [
            sum(i < len(target) and focal[i] == target[i] and focal[i] not in '-X*' for i in flank)
            for flank in flanks
        ]
        pair = dict(
            ortholog=oid,
            state=state,
            map=map_row,
            intron=intron_seqs.get(oid, ''),
            alignment_focal=focal[low:high],
            alignment_target=target[low:high],
            alignment_boundary_column=center - low,
            flank_matches_recomputed=matches,
            flank_available=[len(flank) for flank in flanks],
            correction=corrections.get((tid, oid)),
        )

        # the miniprot hit behind the call, with its predicted introns near the junction
        work_dir = ORTHOLOGY / 'primary/orthology' / genome_id / '_work'
        gffs = list(work_dir.glob(f'{oid}.*.miniprot.gff'))
        if len(gffs) == 1:
            gff = gffs[0]
            if gff not in parsed_gffs:
                parsed_gffs[gff] = parse_miniprot(gff)
                source_hashes[str(gff)] = sha(gff)
            hit = parsed_gffs[gff].get(protein)
            if hit:
                nearby, longest = nearby_prediction(hit, residue)
                pair.update(
                    miniprot_hit=hit,
                    nearby_predicted_introns=nearby,
                    selected_predicted_intron=longest,
                    miniprot_source=str(gff),
                )
        pairs.append(pair)
    return protein, residue, intron_seqs, pairs


def main():
    check_scratch()
    OUT.mkdir(parents=True, exist_ok=True)

    # inputs: genome metadata, primary catalogue, closest comparators and RNA support
    meta = {
        row['genome_id']: row for row in json.loads((ORTHOLOGY / 'inputs/genomes.json').read_text())
    }
    catalogue = {row['locus_id']: row for row in read_tsv(ORTHOLOGY / 'inputs/primary.tsv')}
    rna = {row['locus_id']: row for row in read_tsv(RNA_LOCI) if row['assayed'] == '1'}
    closest_by_genome = defaultdict(list)
    for row in primary_closest():
        closest_by_genome[row['genome_id']].append(row)

    genomes = []
    loci = []
    selection = []
    source_hashes = {}

    # per genome: reviewed states, recall map rows and absence-correction decisions
    for genome_id, closest_rows in closest_by_genome.items():
        reviewed = ORTHOLOGY / 'reviewed/primary/orthology' / genome_id / 'states.tsv'
        map_path = ORTHOLOGY / 'recall/primary/orthology' / genome_id / 'ortholog_map.tsv'
        decisions = ORTHOLOGY / 'corrected/primary/orthology' / genome_id / 'decisions.json'
        for path in (reviewed, map_path, decisions):
            source_hashes[str(path)] = sha(path)
        states_by_locus = defaultdict(list)
        for row in read_tsv(reviewed):
            states_by_locus[row['telotron']].append(row)
        map_rows = defaultdict(list)
        maps = read_tsv(map_path)
        for row in maps:
            map_rows[row['telotron']].append(row)
        corrections = {
            (decision['locus'], decision['ortholog']): decision
            for decision in json.loads(decisions.read_text())['decisions']
        }

        # zero: the genome has no closest NO_INTRON locus
        zero = not any(row['closest_state'] == 'NO_INTRON' for row in closest_rows)
        chosen = [row for row in closest_rows if row['closest_state'] == 'ORDINARY']
        resolved = [
            row
            for row in closest_rows
            if any(s['state'] in RESOLVED_STATES for s in states_by_locus[row['locus_id']])
        ]
        genome_states = [s for row in closest_rows for s in states_by_locus[row['locus_id']]]
        genomes.append(
            dict(
                genome_id=genome_id,
                organism=meta[genome_id]['organism'],
                zero_percent=zero,
                primary=len(closest_rows),
                position_resolved=len(resolved),
                closest_states=dict(Counter(row['closest_state'] for row in closest_rows)),
                comparator_states=dict(Counter(s['state'] for s in genome_states)),
                raw_calls=dict(Counter(row['call'] for row in maps)),
                reviewed_loci=len(chosen),
            )
        )

        parsed_gffs = {}
        for closest in chosen:
            tid = closest['locus_id']
            states = states_by_locus[tid]
            protein, residue, intron_seqs, pairs = locus_pairs(
                tid, genome_id, states, map_rows[tid], corrections, parsed_gffs, source_hashes
            )

            # the focal intron must match the catalogue hash
            focal_intron = intron_seqs.get(genome_id)
            focal_sha = (
                hashlib.sha256(focal_intron.upper().encode()).hexdigest() if focal_intron else ''
            )
            assert focal_sha == catalogue[tid]['intron_sha256'], tid

            loci.append(
                dict(
                    locus_id=tid,
                    genome_id=genome_id,
                    organism=meta[genome_id]['organism'],
                    zero_percent_genome=zero,
                    selection_reason=['all current closest-ordinary loci'],
                    catalogue=catalogue[tid],
                    closest=closest,
                    all_accepted_states=states,
                    all_raw_mapping_calls=map_rows[tid],
                    focal_intron=focal_intron,
                    protein_id=protein,
                    intron_residue=residue,
                    intron_phase=map_rows[tid][0]['intron_phase'],
                    rna=rna.get(tid),
                    pairs=pairs,
                )
            )
            selection.append(
                dict(
                    locus_id=tid,
                    genome_id=genome_id,
                    reason='all current closest-ordinary loci',
                    closest_state=closest['closest_state'],
                )
            )
        print(genome_id, len(chosen), 'loci extracted', flush=True)

    n_pairs = sum(len(locus['pairs']) for locus in loci)
    write_json_gz('evidence.json.gz', dict(genomes=genomes, loci=loci, source_hashes=source_hashes))
    (OUT / 'selection.json').write_text(json.dumps(selection, indent=2) + '\n')
    zero_resolved = sum(
        locus['zero_percent_genome']
        and any(s['state'] in RESOLVED_STATES for s in locus['all_accepted_states'])
        for locus in loci
    )
    complete = dict(
        complete=True,
        genomes=len(genomes),
        loci=len(loci),
        pairs=n_pairs,
        zero_genome_resolved=zero_resolved,
        script_sha256=sha(__file__),
        evidence_sha256=sha(OUT / 'evidence.json.gz'),
        selection_sha256=sha(OUT / 'selection.json'),
    )
    (OUT / 'complete.json').write_text(json.dumps(complete, indent=2) + '\n')
    print('COMPLETE', len(loci), 'loci', flush=True)


if __name__ == '__main__':
    main()
