"""Fetch the DNA of hand-picked alternative comparator introns.

Input config/alternative_requests.json is written by hand: a list of objects,
one per alternative intron, each with
  key        '<locus_id>|<comparator genome_id>' (free text, carried through)
  locus_id   the focal telotron locus the alternative belongs to
  genome_id  the comparator genome holding the alternative intron
  seqid, strand, start, end   its 1-based inclusive coordinates and strand
  motif      the focal lineage telomere motif used to score it (measure_alternatives.py)
Other fields are carried through. The sequence is read from the genome FASTA under data/raw
(via the stage 5 reference cache when one exists).
Output: work/06_exact_position/alternative_sequences.json.
"""

import hashlib
import json

from exact_position import OUT, REPO, GenomeSequence, check_scratch, genome_files


def main():
    check_scratch()
    requests = json.loads((REPO / 'config/alternative_requests.json').read_text())
    for request in requests:
        fasta, _ = genome_files(request['genome_id'])
        assert fasta, request['genome_id']
        genome = GenomeSequence(fasta, {request['seqid']})
        start, end = request['start'] - 1, request['end']
        seq = genome.fetch_stranded(request['seqid'], start, end, request['strand'])
        genome.close()
        assert len(seq) == end - start, request['key']
        request['reference'] = genome.provenance
        request['sequence'] = seq
        request['sequence_sha256'] = hashlib.sha256(seq.encode()).hexdigest()

    (OUT / 'alternative_sequences.json').write_text(json.dumps(requests, indent=2) + '\n')
    print('COMPLETE', len(requests), flush=True)


if __name__ == '__main__':
    main()
