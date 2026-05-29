#!/usr/bin/env python3
"""BLAST telotron linker sequences against genomes.

Stages:
    1. extract per-genome linker FASTAs from architecture TSV
    2. build a per-genome blastn db for every positive species
    3. build a combined "all genomes" blastn db (concatenation)
    4. blast each linker set against its OWN genome  -> hits_own.tsv
    5. blast all linkers against the combined db     -> hits_all.tsv

Linker query header format:
    >{gid}|{seqid}|{start}-{end}|strand|arch|llen={n}

Output TSV columns (outfmt 6 + qcovs + qlen, plus query_genome_id parsed from qseqid):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen query_genome_id

Decisions:
    - dustmasker is OFF (-dust no): linker is non-telomeric by definition but is short,
      and the default dust low-complexity filter eats short repeats we want to see.
    - word_size 11, e-value 1e-5: lenient enough to catch ~70% identity over 30 bp
      (typical interspersed-repeat fragment) but rejects random matches.
    - max_target_seqs 50: for "is this linker a recurring genomic sequence?" we
      want to see multiple hits if they exist.
"""
import argparse
import gzip
import os
import shutil
import subprocess as sp
import sys

import pandas as pd

from _common import slug as _slug, find_genome_fasta


def decompress(src, dst):
    """Decompress .gz to dst if needed; else copy."""
    if src.endswith(".gz"):
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(src, dst)


def run(cmd, **kw):
    print("$", " ".join(cmd), flush=True)
    return sp.run(cmd, check=True, **kw)


BLAST_FIELDS = ("qseqid sseqid pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore qcovs qlen")
HEADER = BLAST_FIELDS.replace(" ", "\t") + "\tquery_genome_id\n"


def parse_qgid(qseqid):
    return qseqid.split("|", 1)[0]


def write_linker_fastas(arch_tsv, outdir, min_linker_len=15):
    """Per-genome linker FASTA. Returns dict gid -> path."""
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(arch_tsv, sep="\t")
    df = df[df.architecture.fillna("").str.contains("linker", na=False)]
    df = df[df.linker_seq.fillna("").str.len() >= min_linker_len]
    out = {}
    for gid, sub in df.groupby("genome_id", sort=False):
        path = f"{outdir}/{_slug(gid)}.fa"
        with open(path, "w") as fh:
            for r in sub.itertuples(index=False):
                header = (f">{gid}|{r.seqid}|{r.start}-{r.end}|{r.strand}"
                          f"|{r.architecture}|llen={r.linker_len}")
                fh.write(header + "\n" + r.linker_seq + "\n")
        out[gid] = path
        print(f"  linker fasta: {gid}  ({len(sub)} linkers)")
    return out


def build_db(fa_path, db_prefix):
    """makeblastdb on a (possibly gz) FASTA; decompresses if needed."""
    if not os.path.exists(db_prefix + ".nhr") and not os.path.exists(db_prefix + ".nsq"):
        if fa_path.endswith(".gz"):
            plain = db_prefix + ".fna"
            decompress(fa_path, plain)
            fa_used = plain
        else:
            fa_used = fa_path
        run(["makeblastdb", "-in", fa_used, "-dbtype", "nucl",
             "-out", db_prefix])


def build_combined_db(linker_gids, refseq_dir, tara_dir, combined_fa, db_prefix):
    """Concatenate every positive-species FASTA (decompressed) into one file, makeblastdb."""
    if os.path.exists(db_prefix + ".nhr") or os.path.exists(db_prefix + ".00.nhr"):
        print(f"combined db already exists at {db_prefix}, skipping")
        return
    with open(combined_fa, "wb") as out_fh:
        for gid in linker_gids:
            try:
                fa = find_genome_fasta(gid, refseq_dir, tara_dir)
            except FileNotFoundError:
                print(f"  skip combined: {gid} fasta not found")
                continue
            opener = gzip.open if fa.endswith(".gz") else open
            with opener(fa, "rb") as fh:
                # prefix contig names with genome_id so we can attribute hits
                for line in fh:
                    if line.startswith(b">"):
                        out_fh.write(b">" + gid.encode() + b"|" + line[1:])
                    else:
                        out_fh.write(line)
    run(["makeblastdb", "-in", combined_fa, "-dbtype", "nucl",
         "-out", db_prefix])


def blastn(query_fa, db_prefix, out_tsv, threads):
    cmd = ["blastn",
           "-query", query_fa,
           "-db", db_prefix,
           "-out", out_tsv,
           "-outfmt", f"6 {BLAST_FIELDS}",
           "-evalue", "1e-5",
           "-word_size", "11",
           "-dust", "no",
           "-max_target_seqs", "50",
           "-num_threads", str(threads)]
    run(cmd)


def annotate_qgid(in_tsv, out_tsv):
    """Prepend column header and append query_genome_id column."""
    with open(in_tsv) as fin, open(out_tsv, "w") as fout:
        fout.write(HEADER)
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            qseqid = line.split("\t", 1)[0]
            fout.write(line + "\t" + parse_qgid(qseqid) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch-tsv", required=True, help="final_telotron_set_architecture.tsv")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--workdir", required=True, help="where to put dbs + queries + raw blast outputs")
    ap.add_argument("--out-own", required=True, help="combined own-genome hits TSV")
    ap.add_argument("--out-all", required=True, help="combined whole-db hits TSV")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-linker-len", type=int, default=15)
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    query_dir = f"{args.workdir}/linker_queries"
    db_dir = f"{args.workdir}/dbs"
    raw_dir = f"{args.workdir}/raw_hits"
    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)

    print("=== 1. linker FASTAs ===")
    gid_to_fa = write_linker_fastas(args.arch_tsv, query_dir, args.min_linker_len)
    if not gid_to_fa:
        print("no linkers; nothing to blast")
        return

    print("\n=== 2. per-genome blastdbs ===")
    gid_to_db = {}
    for gid in gid_to_fa:
        try:
            fa = find_genome_fasta(gid, args.refseq_dir, args.tara_dir)
        except FileNotFoundError:
            print(f"  skip own db: {gid}")
            continue
        db_prefix = f"{db_dir}/{_slug(gid)}"
        build_db(fa, db_prefix)
        gid_to_db[gid] = db_prefix

    print("\n=== 3. combined db ===")
    combined_fa = f"{db_dir}/all_genomes.fna"
    combined_db = f"{db_dir}/all_genomes"
    build_combined_db(list(gid_to_fa.keys()), args.refseq_dir, args.tara_dir,
                      combined_fa, combined_db)
    # the concatenated fasta can be huge; drop it once db is built
    if os.path.exists(combined_fa) and os.path.exists(combined_db + ".nhr"):
        os.remove(combined_fa)

    print("\n=== 4. blast linkers vs own genome ===")
    own_hits = []
    for gid, qfa in gid_to_fa.items():
        if gid not in gid_to_db:
            continue
        raw = f"{raw_dir}/{_slug(gid)}_own.tsv"
        blastn(qfa, gid_to_db[gid], raw, args.threads)
        own_hits.append(raw)

    with open(args.out_own, "w") as out_fh:
        out_fh.write(HEADER)
        for raw in own_hits:
            with open(raw) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line: continue
                    out_fh.write(line + "\t" + parse_qgid(line.split("\t", 1)[0]) + "\n")
    print(f"wrote {args.out_own}")

    print("\n=== 5. blast all linkers vs combined db ===")
    all_query = f"{query_dir}/_all_linkers.fa"
    with open(all_query, "wb") as out_fh:
        for qfa in gid_to_fa.values():
            with open(qfa, "rb") as fh:
                shutil.copyfileobj(fh, out_fh)
    raw_all = f"{raw_dir}/all_vs_combined.tsv"
    blastn(all_query, combined_db, raw_all, args.threads)
    annotate_qgid(raw_all, args.out_all)
    print(f"wrote {args.out_all}")


if __name__ == "__main__":
    sys.exit(main())
