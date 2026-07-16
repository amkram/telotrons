#!/usr/bin/env python3
"""Deduplicate doubly-assembled exons and close gene copies.

Rationale
    Telotron-scanning a single species sometimes returns the same locus more
    than once: (a) the assembly has duplicated contigs (haplotype splits,
    misassembled tandem regions); (b) the locus sits in a paralog whose
    flanking exons are nearly identical to another paralog. Without
    deduplication, downstream alignments / architecture stats double-count the
    locus and inflate apparent counts.

Method
    For every species in the input TSV:
      1. Extract 250 bp upstream and 250 bp downstream flanks per locus
         (spliced/display orientation) and write each as a separate FASTA record.
      2. blastn all-vs-all of the flanks (own DB, -dust no, word_size 11).
      3. Two loci A and B are duplicates if BOTH their upstream flanks AND
         their downstream flanks share an alignment of >= MIN_LEN bp at
         >= MIN_PIDENT identity (Ling 2007 — Apicomplexan multigene families
         have ≥95% identity in flanks; requiring both sides simultaneously
         prevents merging recent paralogs that independently acquired a
         telotron). Self-hits (A==B) and cross-side hits (A.upstream vs
         B.downstream) are ignored.
      4. Connected components in the duplicate graph collapse to one
         representative — the locus with the longest intron (tie-break: first
         occurrence). All other members are dropped.

Outputs
    --out-final  deduplicated TSV (subset of input rows)
    --out-log    per-locus action log with group id and representative
"""
import argparse
import os
import subprocess as sp
import sys
import tempfile

import pandas as pd

from _common import rc, slug as _slug, load_fasta, find_genome_fasta

MIN_PIDENT = 95.0
MIN_LEN = 100
FLANK = 250


def flanks_for_locus(row, chrom, flank=FLANK):
    """Return (upstream_seq, downstream_seq) in spliced/display orientation."""
    L = len(chrom)
    ls = max(0, row.start - 1 - flank)
    le = row.start - 1
    rs = row.end
    re_ = min(L, row.end + flank)
    if row.strand == "-":
        up = rc(chrom[rs:re_])
        dn = rc(chrom[ls:le])
    else:
        up = chrom[ls:le]
        dn = chrom[rs:re_]
    return up, dn


def union_find_groups(n, edges):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    for a, b in edges:
        union(a, b)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def dedup_one_species(gid, sub, chrom_seqs, tmpdir, threads):
    """Return (keep_mask: list[bool], group_ids: list[int], rep_indices: list[int]).

    sub is the per-species dataframe slice. Indexing is 0..len(sub)-1.
    """
    n = len(sub)
    if n <= 1:
        return [True] * n, [0] * n, list(range(n))

    fa_path = f"{tmpdir}/{_slug(gid)}_flanks.fa"
    db_prefix = f"{tmpdir}/{_slug(gid)}_db"
    out_blast = f"{tmpdir}/{_slug(gid)}_blast.tsv"

    # Write upstream + downstream as separate records: id = "LOCUSidx|U" or "|D"
    with open(fa_path, "w") as fh:
        for i, r in enumerate(sub.itertuples(index=False)):
            chrom = chrom_seqs.get(r.seqid, "")
            if not chrom:
                continue
            up, dn = flanks_for_locus(r, chrom)
            if up:
                fh.write(f">L{i}|U\n{up}\n")
            if dn:
                fh.write(f">L{i}|D\n{dn}\n")

    if os.path.getsize(fa_path) == 0:
        return [True] * n, list(range(n)), list(range(n))

    sp.run(["makeblastdb", "-in", fa_path, "-dbtype", "nucl", "-out", db_prefix],
           check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    # Drop -max_target_seqs: for species with several thousand loci
    # (E. tenella ~3000; TARA_PSW_86_MAG_00284 ~3100) a per-query hit list
    # exceeding the cap silently loses the true duplicate partner AND
    # biases which hits win (BLAST's cap is order-of-encounter, not by
    # e-value). We downstream-filter by pident/len anyway.
    sp.run(["blastn", "-query", fa_path, "-db", db_prefix, "-out", out_blast,
            "-outfmt", "6 qseqid sseqid pident length",
            "-evalue", "1e-10", "-word_size", "11", "-dust", "no",
            "-num_threads", str(threads)],
           check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

    # Collect per-side duplicate pairs separately; an edge is admitted only when
    # BOTH the upstream and downstream flanks of the same locus pair are duplicates
    # (Apicomplexan SAG/multigene families often share ≥95% flank identity on one
    # side without being the same locus — requiring both sides preserves recent
    # paralog insertions; Ling 2007).
    pairs_U = set()
    pairs_D = set()
    with open(out_blast) as fh:
        for line in fh:
            qid, sid, pid, ln = line.rstrip("\n").split("\t")
            qi_s, qside = qid.split("|")
            si_s, sside = sid.split("|")
            if qside != sside:
                continue                # only same-side hits count
            try:
                qi, si = int(qi_s[1:]), int(si_s[1:])
            except ValueError:
                continue                # malformed blast record id
            if qi == si:
                continue                # self
            if float(pid) < MIN_PIDENT or int(ln) < MIN_LEN:
                continue
            key = (qi, si) if qi < si else (si, qi)
            if qside == "U":
                pairs_U.add(key)
            else:
                pairs_D.add(key)
    edges = [(a, b) for (a, b) in pairs_U & pairs_D]
    # One-side-only edge diagnostic: pairs where a single flank shares
    # >=MIN_LEN/MIN_PIDENT but not both. These are potential under-collapse
    # candidates (segmental-duplication artifacts collapsed on one side only).
    one_side_only = len((pairs_U | pairs_D) - (pairs_U & pairs_D))
    print(f"    [dedup] both-flank edges: {len(edges)}; one-side-only (potential under-collapse): {one_side_only}",
          flush=True)

    groups = union_find_groups(n, edges)
    keep = [False] * n
    group_id = [0] * n
    rep_of = [0] * n
    for gid_idx, members in enumerate(groups):
        # representative: longest intron; tie-break by lowest index
        rep = max(members, key=lambda i: (sub.iloc[i].intron_len, -i))
        for m in members:
            group_id[m] = gid_idx
            rep_of[m] = rep
        keep[rep] = True
    return keep, group_id, rep_of


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True, help="final_telotron_set.tsv (post-filter)")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--out-final", required=True, help="deduplicated TSV")
    ap.add_argument("--out-log", required=True, help="per-locus dedup log TSV")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    df = pd.read_csv(args.final, sep="\t")
    keep_mask_all = []
    log_rows = []

    with tempfile.TemporaryDirectory(prefix="dedup_telotrons_") as tmpdir:
        for gid, sub in df.groupby("genome_id", sort=False):
            sub = sub.reset_index(drop=False).rename(columns={"index": "_orig_idx"})
            try:
                fa = find_genome_fasta(gid, args.refseq_dir, args.tara_dir)
                seqs = load_fasta(fa)
            except FileNotFoundError:
                print(f"skip dedup for {gid} (no FASTA); keeping all {len(sub)} loci")
                for i, r in sub.iterrows():
                    log_rows.append({"genome_id": gid, "_orig_idx": r._orig_idx,
                                     "seqid": r.seqid, "start": r.start, "end": r.end,
                                     "group_id": i, "representative_idx": r._orig_idx,
                                     "kept": True, "reason": "no_fasta"})
                keep_mask_all.extend([True] * len(sub))
                continue

            keep, group_id, rep_of = dedup_one_species(gid, sub, seqs, tmpdir, args.threads)
            n_in = len(sub)
            n_out = sum(keep)
            print(f"  {gid}: {n_in} -> {n_out} loci  ({n_in - n_out} dropped as duplicates)")

            for i, r in sub.iterrows():
                rep_orig = sub.iloc[rep_of[i]]._orig_idx
                log_rows.append({"genome_id": gid, "_orig_idx": r._orig_idx,
                                 "seqid": r.seqid, "start": r.start, "end": r.end,
                                 "group_id": group_id[i],
                                 "representative_idx": rep_orig,
                                 "kept": keep[i],
                                 "reason": "kept" if keep[i] else "duplicate"})
            keep_mask_all.extend(keep)

    log_df = pd.DataFrame(log_rows)
    log_df.to_csv(args.out_log, sep="\t", index=False)

    # Apply keep_mask_all in the order rows were processed
    dedup_df = df.iloc[[r["_orig_idx"] for r in log_rows if r["kept"]]].copy()
    dedup_df.to_csv(args.out_final, sep="\t", index=False)
    print(f"\ndedup input:  {len(df)} loci")
    print(f"dedup output: {len(dedup_df)} loci  ({len(df) - len(dedup_df)} dropped)")
    print(f"wrote {args.out_final}")
    print(f"wrote {args.out_log}")


if __name__ == "__main__":
    sys.exit(main())
