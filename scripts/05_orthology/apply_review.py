"""Apply the curated pair review (config/orthology_review.tsv) to the corrected states of every
genome in both tiers, writing reviewed/<stage>/orthology/<genome>/states.tsv.

Each review row names a telotron x comparator pair, the state it must have after the absence
correction, the state that replaces it and why: UNRESOLVED_MAPPING for comparator hits shown to
join a remote exon fragment to the conserved core, AMBIGUOUS for a comparator intron that is
itself an unlisted telomeric repeat variant. UNRESOLVED_MAPPING also replaces the aligner call.
"""

from collections import Counter

from common import REPO, STAGES, genomes, header, orthology_dir, read, write

REVIEW = REPO / 'config/orthology_review.tsv'


def main():
    review = read(REVIEW)
    changes = {(row['stage'], row['telotron'], row['ortholog']): row for row in review}
    assert len(changes) == len(review), 'duplicate review pairs'

    applied = set()
    for stage in STAGES:
        for meta in genomes():
            corrected = orthology_dir(stage, meta['genome_id'], 'corrected') / 'states.tsv'
            rows = read(corrected)
            for row in rows:
                key = (stage, row['telotron'], row['ortholog'])
                if key not in changes:
                    continue
                change = changes[key]
                assert row['state'] == change['original_state'], (key, row['state'])
                row['state'] = change['replacement_state']
                if change['replacement_state'] == 'UNRESOLVED_MAPPING':
                    row['aligner_call'] = 'UNRESOLVED_MAPPING'
                applied.add(key)
            reviewed = orthology_dir(stage, meta['genome_id'], 'reviewed') / 'states.tsv'
            write(reviewed, rows, header(corrected))

    # every review row must match exactly one pair
    assert applied == set(changes), sorted(set(changes) - applied)
    counts = Counter(row['replacement_state'] for row in review)
    print('REVIEWED', len(applied), 'pairs', dict(counts), flush=True)


if __name__ == '__main__':
    main()
