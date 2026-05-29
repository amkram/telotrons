#!/usr/bin/env python3
"""BLAST the telomerase reference sequences (telomerase.us) against each genome.

For every (query-set, genome) pair, run BLAST in permissive mode with `outfmt 6`
and interleave each tabular row with the ASCII pairwise alignment of that HSP.

Query → BLAST program  → output filename
  tr.nucl     blastn   → {genome}/tr_nucl.tsv
  tert.nucl   blastn   → {genome}/tert_nucl.tsv
  tert.prot   tblastn  → {genome}/tert_prot.tsv
  other.nucl  blastn   → {genome}/other_nucl.tsv
  other.prot  tblastn  → {genome}/other_prot.tsv

Parameters (permissive vs defaults, but tight enough to finish quickly):
  blastn:   -evalue 1e-3 -word_size 11
  tblastn:  -evalue 1e-3 -word_size 5
DUST and SEG masking left at BLAST defaults (enabled).

Per HSP we emit:
  <tabular line: outfmt 6 standard cols + qseq + sseq>
  qseq:  <aligned query>
  match: <midline (| for matches, ' ' otherwise)>
  sseq:  <aligned subject>
  (blank line)
"""
import argparse
import gzip
import os
import shutil
import subprocess as sp
import sys
import tempfile

from _common import find_genome_fasta

QUERIES = [
    # (label, fasta, program, dbtype)
    ("tr_nucl",    "tr.nucl.fa",    "blastn",  "nucl"),
    ("tert_nucl",  "tert.nucl.fa",  "blastn",  "nucl"),
    ("tert_prot",  "tert.prot.fa",  "tblastn", "prot"),  # prot vs nucl db (translated)
    ("other_nucl", "other.nucl.fa", "blastn",  "nucl"),
    ("other_prot", "other.prot.fa", "tblastn", "prot"),
]

# outfmt 6 columns we request (12 standard + qseq + sseq). The qseq/sseq trail
# the standard 12-col schema so external tools still see the usual layout.
OUTFMT = ("6 qseqid sseqid pident length mismatch gapopen "
          "qstart qend sstart send evalue bitscore qseq sseq")
HEADER = ("# qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\t"
          "qstart\tqend\tsstart\tsend\tevalue\tbitscore\tqseq\tsseq")


def decompress_if_needed(src, tmpdir):
    if src.endswith(".gz"):
        out = os.path.join(tmpdir, os.path.basename(src)[:-3])
        with gzip.open(src, "rb") as fin, open(out, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        return out
    return src


def build_blastdb(fa, dbtype, tmpdir):
    """Make a temporary blastdb under tmpdir; return its prefix."""
    prefix = os.path.join(tmpdir, "genome")
    sp.run(["makeblastdb", "-in", fa, "-dbtype", dbtype, "-out", prefix,
            "-parse_seqids"],
           check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    return prefix


def midline(qseq, sseq, is_protein):
    """Match line: '|' for identity, '+' for positive (protein), ' ' otherwise."""
    out = []
    for q, s in zip(qseq, sseq):
        if q == "-" or s == "-":
            out.append(" ")
        elif q.upper() == s.upper():
            out.append("|")
        else:
            out.append(" ")
    return "".join(out)


def write_block(out_fh, row, is_protein):
    """Tabular row + 3-line alignment + blank separator."""
    cols = row.rstrip("\n").split("\t")
    out_fh.write("\t".join(cols[:12]) + "\n")  # outfmt 6 std (drop qseq/sseq from row)
    if len(cols) >= 14:
        qseq, sseq = cols[12], cols[13]
        out_fh.write(f"  qseq:  {qseq}\n")
        out_fh.write(f"  match: {midline(qseq, sseq, is_protein)}\n")
        out_fh.write(f"  sseq:  {sseq}\n")
    out_fh.write("\n")


def run_blast(prog, query_fa, db_prefix, threads, is_protein):
    """Stream HSPs in our extended outfmt 6. Returns full stdout text."""
    if prog == "blastn":
        params = ["-evalue", "1e-3", "-word_size", "11"]
    else:  # tblastn
        params = ["-evalue", "1e-3", "-word_size", "5"]
    cmd = [prog, "-query", query_fa, "-db", db_prefix,
           "-outfmt", OUTFMT, "-num_threads", str(threads)] + params
    res = sp.run(cmd, check=True, capture_output=True, text=True)
    return res.stdout


def process_one_genome(gid, fa_src, db_dir, outdir, queries, threads):
    """Build a blastdb for one genome and run every query set against it."""
    with tempfile.TemporaryDirectory(prefix=f"blast_{gid}_") as tmpdir:
        fa = decompress_if_needed(fa_src, tmpdir)
        db_prefix = build_blastdb(fa, "nucl", tmpdir)
        os.makedirs(outdir, exist_ok=True)
        for label, q_fa_name, prog, _ in queries:
            q_path = os.path.join(db_dir, q_fa_name)
            if not (os.path.exists(q_path) and os.path.getsize(q_path) > 0):
                print(f"  skip {label}: {q_path} missing/empty", flush=True)
                continue
            out_path = os.path.join(outdir, f"{label}.tsv")
            txt = run_blast(prog, q_path, db_prefix, threads,
                            is_protein=(prog == "tblastn"))
            with open(out_path, "w") as out:
                out.write(HEADER + "\n")
                n = 0
                for line in txt.splitlines():
                    if not line.strip():
                        continue
                    write_block(out, line, is_protein=(prog == "tblastn"))
                    n += 1
            print(f"  {label}: {n} HSPs → {out_path}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-id", required=True,
                    help="single genome to BLAST against (one per script call)")
    ap.add_argument("--db-dir", required=True,
                    help="dir with tr.nucl.fa, tert.nucl.fa, tert.prot.fa, ...")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True,
                    help="per-genome output dir (created if needed)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    gid = args.genome_id
    fa_src = find_genome_fasta(gid, args.refseq_dir, args.tara_dir, required=False)
    if not fa_src:
        print(f"skip {gid}: FASTA not found", flush=True)
        sys.exit(0)
    print(f"[{gid}]  {fa_src}", flush=True)
    process_one_genome(gid, fa_src, args.db_dir, args.outdir,
                       QUERIES, threads=args.threads)


if __name__ == "__main__":
    sys.exit(main())
