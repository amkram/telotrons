"""Re-call one genome's comparisons at >=12 matches per 20-residue flank (the reported gate)."""

import json
from collections import defaultdict

from classify_states import CALLER, classify
from common import MIN_FLANK_MATCHES, RECALL_POLICY, REPO, WORK, dump, orthology_dir, parse_stage
from common import read, sha
from common import task_genome, use_scratch, write
from indexed_reference import Reference
from telotron_ortholog_align import (
    FLANK_AA,
    INTRON_TOL,
    MAX_ALT_LOCI,
    MIN_IDENT,
    intron_dna,
    introns_elsewhere,
    longest_intron,
    nearby_introns,
    parse_miniprot,
    residue_position,
    site_clear,
    structurally_absent,
    telomere_frac,
)
from _common import find_genome_files, locus_stub, read_aln, read_fasta


def add_flank_matches(rows, proteins):
    """Add flank residue counts to map rows. Returns rows >= MIN_IDENT by comparator."""
    evaluated = defaultdict(list)
    for row in rows:
        residue = int(row['intron_residue'])
        protein_len = len(proteins[row['protein']])
        assert 0 <= residue <= protein_len
        available5 = min(FLANK_AA, residue)
        available3 = min(FLANK_AA, protein_len - residue)

        # the aligner measured flank mismatches only above its identity floor
        above_floor = float(row['pct_identity']) >= MIN_IDENT
        row.update(
            flank5_available=available5,
            flank3_available=available3,
            flank5_actual_matches=available5 - int(row['flank5_mm']) if above_floor else '',
            flank3_actual_matches=available3 - int(row['flank3_mm']) if above_floor else '',
            orthology_policy=RECALL_POLICY,
        )
        if above_floor:
            evaluated[row['ortholog']].append(row)
    return evaluated


def recall_pair(row, hit, contig, motif):
    """Re-call one comparison under the 12-match gate and return its FILL intron DNA, or ''."""
    residue = int(row['intron_residue'])
    gate = (
        row['flank5_actual_matches'] >= MIN_FLANK_MATCHES
        and row['flank3_actual_matches'] >= MIN_FLANK_MATCHES
        and row['rbh_pass'] == '1'
        and int(row['n_alternative_loci']) <= MAX_ALT_LOCI
    )

    # same structural rules as the aligner
    candidates = nearby_introns(hit, residue)
    struct_absent = False
    absent = False
    seq = ''
    if candidates and gate:
        seq = intron_dna(contig, hit, *longest_intron(candidates))
        assert seq
    elif not candidates:
        struct_absent = structurally_absent(hit, residue)
        if gate and struct_absent and introns_elsewhere(hit, residue):
            position = residue_position(hit, residue)
            absent = bool(contig) and position is not None and site_clear(contig, position)

    # precedence as in the aligner: FILL, ABSENT, POOR_FLANK, UNCERTAIN, NOCOV
    if seq:
        row['call'] = 'FILL'
    elif absent:
        row['call'] = 'ABSENT'
    elif not gate and (candidates or struct_absent):
        row['call'] = 'POOR_FLANK'
    elif hit['qmin'] - INTRON_TOL <= residue <= hit['qmax'] + INTRON_TOL:
        row['call'] = 'UNCERTAIN'
    else:
        row['call'] = 'NOCOV'
    row['ortho_intron_len'] = len(seq)
    row['ortho_intron_telofrac'] = f"{telomere_frac(seq, motif):.3f}" if seq else '0.000'
    return seq


def check_against_baseline(rows, baseline, states):
    """Kept states pass the gate and match the baseline. Lost baseline states had short flanks."""
    by_pair = {(row['telotron'], row['ortholog']): row for row in rows}
    kept = {(state['telotron'], state['ortholog']): state for state in states}
    for pair, state in kept.items():
        row = by_pair[pair]
        matches = min(row['flank5_actual_matches'], row['flank3_actual_matches'])
        assert matches >= MIN_FLANK_MATCHES
        assert row['rbh_pass'] == '1'
        assert int(row['n_alternative_loci']) <= MAX_ALT_LOCI
        assert row['frameshift_at_junction'] == '0'
        if pair in baseline:
            assert state['state'] == baseline[pair]['state'], ('Changed state', pair)

    removed = baseline.keys() - kept.keys()
    for pair in removed:
        row = by_pair[pair]
        matches = min(row['flank5_actual_matches'], row['flank3_actual_matches'])
        assert matches < MIN_FLANK_MATCHES, ('Unexpected baseline loss', pair, row)
    return dict(
        baseline_states=len(baseline),
        baseline_retained=len(baseline.keys() & kept.keys()),
        new_states=len(kept.keys() - baseline.keys()),
        removed_short_flanks=len(removed),
    )


def main():
    use_scratch()
    stage = parse_stage(__doc__)
    meta = task_genome()
    genome_id = meta['genome_id']
    source_dir = orthology_dir(stage, genome_id)
    out_dir = orthology_dir(stage, genome_id, 'recall')
    out_dir.mkdir(parents=True, exist_ok=True)

    # the source run must be complete, classified and unchanged
    certificate = json.loads((source_dir / 'complete.json').read_text())
    assert certificate['complete']
    assert certificate.get('composition'), ('comparator states not classified', genome_id)
    assert sha(source_dir / 'ortholog_map.tsv') == certificate['map_sha256']
    assert sha(source_dir / 'states.tsv') == certificate['states_sha256']
    certificate_sha = sha(source_dir / 'complete.json')

    # skip a finished recall, which must match this source, script and outputs
    if (out_dir / 'complete.json').exists():
        prior = json.loads((out_dir / 'complete.json').read_text())
        assert prior['source_certificate_sha256'] == certificate_sha
        assert prior['rerun_script_sha256'] == sha(__file__)
        assert prior['states_sha256'] == sha(out_dir / 'states.tsv')
        assert prior['map_sha256'] == sha(out_dir / 'ortholog_map.tsv')
        print('ALREADY COMPLETE', stage, genome_id, flush=True)
        return

    # source map rows, baseline states, mapped catalogue and focal proteins
    rows = read(source_dir / 'ortholog_map.tsv')
    states = read(source_dir / 'states.tsv')
    baseline = {(state['telotron'], state['ortholog']): state for state in states}
    mapped_catalogue = WORK / stage / 'host' / genome_id / 'mapped_catalogue.tsv'
    catalogue = {row['locus_id']: row for row in read(mapped_catalogue)}
    assert set(catalogue) == {row['telotron'] for row in rows}
    proteins = read_fasta(str(source_dir / '_work/focal_host_proteins.faa'))
    rows_by_comparator = add_flank_matches(rows, proteins)

    # focal intron DNA from the source alignments (comparator FILL DNA is added below)
    intron_seqs = defaultdict(dict)
    for tid in catalogue:
        path = source_dir / 'intron_dna' / f'{locus_stub(tid)}.aln.fasta'
        focal = {name.split('|')[0]: s.replace('-', '') for name, s in read_aln(path)}
        assert focal.get(genome_id), ('Missing focal sequence', tid)
        intron_seqs[tid][genome_id] = focal[genome_id]

    # re-call each comparator's rows from its miniprot hits
    for oid, comparator_rows in sorted(rows_by_comparator.items()):
        gffs = list((source_dir / '_work').glob(f'{oid}.*.miniprot.gff'))
        assert len(gffs) == 1, (genome_id, oid, gffs)
        hits = parse_miniprot(gffs[0])
        fasta, _ = find_genome_files(
            oid, str(REPO / 'data/raw/refseq'), str(REPO / 'data/raw/tara')
        )
        reference = Reference(fasta, build=False)
        for row in comparator_rows:
            hit = hits[row['protein']]
            contig = reference.get(hit['contig'], '')
            motif = catalogue[row['telotron']]['motif']
            seq = recall_pair(row, hit, contig, motif)
            if seq:
                intron_seqs[row['telotron']][oid] = seq
        print('RECALLED', stage, genome_id, oid, len(comparator_rows), flush=True)

    # raw intron DNA per locus for the state caller (--raw-dna), reusing the source protein MSAs
    (out_dir / 'intron_dna_raw').mkdir(exist_ok=True)
    for tid, seqs in intron_seqs.items():
        fasta_text = ''.join(f'>{oid}|raw_intron|len={len(s)}\n{s}\n' for oid, s in seqs.items())
        (out_dir / 'intron_dna_raw' / f'{locus_stub(tid)}.fasta').write_text(fasta_text)
    if not (out_dir / 'protein_msa').exists():
        (out_dir / 'protein_msa').symlink_to(source_dir / 'protein_msa', target_is_directory=True)
    write(out_dir / 'ortholog_map.tsv', rows)

    # classify, then check every kept state against the new gate and the baseline
    composition = classify(stage, meta, recall=True)
    recalled = read(out_dir / 'states.tsv')
    changes = check_against_baseline(rows, baseline, recalled)

    # the source run must still be unchanged
    assert sha(source_dir / 'complete.json') == certificate_sha
    assert sha(source_dir / 'states.tsv') == certificate['states_sha256']

    dump(
        out_dir / 'complete.json',
        dict(
            complete=True,
            genome_id=genome_id,
            stage=stage,
            policy=RECALL_POLICY,
            minimum_actual_matches_per_flank=MIN_FLANK_MATCHES,
            flank_window_aa=FLANK_AA,
            mapped_loci=len(catalogue),
            comparators=len(meta['comparators']),
            states=len(recalled),
            source_certificate_sha256=certificate_sha,
            rerun_script_sha256=sha(__file__),
            map_sha256=sha(out_dir / 'ortholog_map.tsv'),
            states_sha256=sha(out_dir / 'states.tsv'),
            caller_sha256=sha(CALLER),
            composition=composition,
            scratch=use_scratch(),
            **changes,
        ),
    )
    print('ORTHOLOGY COMPLETE', stage, genome_id, len(recalled), changes, flush=True)


if __name__ == '__main__':
    main()
