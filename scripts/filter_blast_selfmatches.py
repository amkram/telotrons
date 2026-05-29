#!/usr/bin/env python3
"""Filter self-matches from telotron linker BLAST output.

A hit is a self-match when the subject overlaps the genomic locus the query
came from.  Query IDs have the form:

    genome_id|seqid|start-end|strand|...

Subject IDs are either bare seqid (own-genome runs) or genome_id|seqid
(all-vs-combined runs).  We compare only the seqid component.

Overlap condition (0-based half-open after converting 1-based BLAST coords):
    not (send < locus_start or sstart > locus_end)

Usage:
    python filter_blast_selfmatches.py hits.tsv > filtered.tsv
    cat hits.tsv | python filter_blast_selfmatches.py - > filtered.tsv

Input: tabular BLAST (-outfmt 6), any number of trailing extra columns.
"""
import sys


def parse_query(qid):
    """Return (seqid, locus_start, locus_end) from a query ID, or None."""
    parts = qid.split("|")
    if len(parts) < 3:
        return None
    seqid = parts[1]
    try:
        start, end = parts[2].split("-")
        return seqid, int(start), int(end)
    except ValueError:
        return None


def subject_seqid(sid):
    """Strip leading genome_id| prefix if present (all-vs-combined format)."""
    parts = sid.split("|")
    return parts[-1] if len(parts) >= 2 else parts[0]


def is_self(qid, sid, sstart, send):
    parsed = parse_query(qid)
    if parsed is None:
        return False
    q_seqid, locus_start, locus_end = parsed
    if subject_seqid(sid) != q_seqid:
        return False
    # BLAST coords are 1-based; treat sstart/send as an inclusive interval.
    s0 = min(sstart, send)
    s1 = max(sstart, send)
    return not (s1 < locus_start or s0 > locus_end)


def main():
    src = open(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "-" else sys.stdin
    kept = dropped = 0
    for line in src:
        if not line.strip() or line.startswith("#"):
            sys.stdout.write(line)
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 10:
            sys.stdout.write(line)
            continue
        qid, sid = cols[0], cols[1]
        try:
            sstart, send = int(cols[8]), int(cols[9])
        except ValueError:
            sys.stdout.write(line)
            continue
        if is_self(qid, sid, sstart, send):
            dropped += 1
        else:
            kept += 1
            sys.stdout.write(line)
    print(f"kept {kept}  dropped {dropped} self-matches", file=sys.stderr)


if __name__ == "__main__":
    main()
