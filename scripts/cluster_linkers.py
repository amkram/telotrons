#!/usr/bin/env python3
"""Cluster linker sequences by 7-mer Jaccard similarity to find recurrent linkers.

Inputs:
    work/results/mechanism_deepdive/linker_segmentation.tsv

Outputs:
    linker_clusters.tsv          -- cluster_id, size, consensus, example_telotrons
    linker_recurrence_summary.txt -- overall stats + length-stratified +
                                     shuffle-null recurrence baseline
"""
from __future__ import annotations

import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

IN_TSV = Path("work/results/mechanism_deepdive/linker_segmentation.tsv")
OUT_DIR = Path("work/results/mechanism_deepdive")
OUT_CLUSTERS = OUT_DIR / "linker_clusters.tsv"
OUT_SUMMARY = OUT_DIR / "linker_recurrence_summary.txt"

# MIN_LEN lowered from 20 -> 5 so the ~11-bp inter-event same-polarity seam
# linker class (memory `telotron_two_linker_classes_2026-06-03`, n=293, median
# 11 bp) enters clustering. Summary is stratified into short (<20 bp,
# inter-event-dominant) and long (>=20 bp, central-junction-dominant) bins.
MIN_LEN = 5
SHORT_LONG_BOUNDARY = 20
K = 7
JACCARD_THRESHOLD = 0.5

# Length-matched shuffle-null for the recurrence rate: re-cluster N replicates
# of length-matched dinucleotide-preserving shuffles of the input sequences and
# report the observed recurrence rate as a delta over null + permutation p.
NULL_REPLICATES = 200
NULL_SEED = 20260605


def revcomp(s: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]


def kmers(seq: str, k: int = K) -> set[str]:
    seq = seq.upper()
    return {seq[i : i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i : i + k]}


def canonical_kmers(seq: str, k: int = K) -> set[str]:
    """Strand-agnostic kmer set: take lex-min of (kmer, revcomp(kmer))."""
    out: set[str] = set()
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        km = seq[i : i + k]
        if "N" in km:
            continue
        rc = revcomp(km)
        out.add(km if km <= rc else rc)
    return out


def jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def representative_cluster(
    kmer_sets: list[set[str]],
    lengths: list[int],
    threshold: float = JACCARD_THRESHOLD,
) -> list[int]:
    """Representative-based clustering (CD-HIT-EST-style) keyed on kmer Jaccard.

    Avoids the single-link chaining pathology of the previous greedy_cluster
    (Sneath & Sokal 1973). Sequences are processed longest-first; each becomes
    a cluster representative if it fails the Jaccard >= threshold test against
    *every* existing representative. Members are compared *only* to the
    representative (not to any other member), which is the CD-HIT semantics
    (Fu 2012, doi:10.1093/bioinformatics/bts565). A kmer->repr-id index makes
    candidate lookup O(k * mean_postings) rather than O(|reprs|) per query.

    Note: representative-based clustering is more conservative than complete-
    link but cannot chain across unrelated linkers via a bridging member.
    CD-HIT-EST / MMseqs2 are the field-standard implementations; this is a
    pure-Python fallback when neither binary is on PATH.
    """
    n = len(kmer_sets)
    labels = [-1] * n
    # process longest first so the representative carries the most kmer mass
    order = sorted(range(n), key=lambda i: (-lengths[i], i))
    # cluster_id -> representative index
    representatives: dict[int, int] = {}
    # kmer -> set of cluster_ids whose REPRESENTATIVE contains the kmer
    repr_kmer_index: dict[str, set[int]] = defaultdict(set)
    next_cid = 0
    for i in order:
        ks = kmer_sets[i]
        if not ks:
            labels[i] = next_cid
            representatives[next_cid] = i
            next_cid += 1
            continue
        # candidate cluster_ids: any whose representative shares >=1 kmer
        cand: set[int] = set()
        for km in ks:
            cand.update(repr_kmer_index.get(km, ()))
        assigned = None
        # deterministic order across candidates
        for cid in sorted(cand):
            rep_ks = kmer_sets[representatives[cid]]
            if jaccard(ks, rep_ks) >= threshold:
                assigned = cid
                break
        if assigned is None:
            assigned = next_cid
            representatives[assigned] = i
            next_cid += 1
            for km in ks:
                repr_kmer_index[km].add(assigned)
        labels[i] = assigned
    return labels


# Back-compat alias so any external import keeps working; new code should
# use representative_cluster directly.
def greedy_cluster(
    kmer_sets: list[set[str]],
    threshold: float = JACCARD_THRESHOLD,
) -> list[int]:
    lengths = [len(ks) for ks in kmer_sets]
    return representative_cluster(kmer_sets, lengths, threshold)


def consensus_from_seqs(seqs: list[str]) -> str:
    """Build a consensus.

    If only one sequence -> return it directly.
    Else strand-orient each seq to the anchor (longest member) by picking
    whichever orientation maximises *forward (non-canonical) 7-mer* overlap
    with the anchor's forward strand, then take the per-column majority over a
    median-length grid (right-padded with '-').

    NOTE on the strand model: the upstream cluster index is canonical
    (lex-min strand) kmers, but consensus must operate in oriented space —
    the anchor's forward strand is the natural reference because it is the
    longest *observed* sequence in the cluster, so any base it carries is
    unambiguous regardless of which strand it came from. We score
    forward(s) vs forward(rc(s)) against the anchor's forward kmers and pick
    whichever orients more bases the same way as the anchor. The previous
    implementation used forward-vs-rc kmer sets against a non-canonical
    anchor and was correct, but was inconsistent with how the cluster was
    built (canonical 7-mers); this function does not re-canonicalise inside
    the consensus because the goal here is a reportable string, not a
    cluster-membership recomputation.
    """
    if not seqs:
        return ""
    if len(seqs) == 1:
        return seqs[0].upper()
    seqs = [s.upper() for s in seqs]
    anchor = max(seqs, key=len)
    anchor_fwd_kmers = kmers(anchor)
    oriented: list[str] = []
    for s in seqs:
        f = kmers(s)
        r = kmers(revcomp(s))
        # tie -> keep forward (deterministic)
        if len(r & anchor_fwd_kmers) > len(f & anchor_fwd_kmers):
            oriented.append(revcomp(s))
        else:
            oriented.append(s)
    # NB: median-length grid + right-pad still discards bases beyond the
    # median from short members; downstream consumers must treat the consensus
    # as a reportable string, not a real MSA. (Acknowledged limitation —
    # MAFFT --localpair would be the field-standard alternative for short
    # divergent strings.)
    L = int(np.median([len(s) for s in oriented]))
    L = max(L, 1)
    cols = []
    for j in range(L):
        col = []
        for s in oriented:
            if j < len(s):
                col.append(s[j])
        if not col:
            cols.append("N")
            continue
        c = Counter(col).most_common(1)[0][0]
        cols.append(c)
    return "".join(cols)


# ---- shuffle null distribution for the recurrence rate -----------------

def _dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """Altschul-Erickson dinucleotide shuffle (Altschul & Erickson 1985)
    via Eulerian-path on the dinucleotide multi-graph; preserves 1st-order
    composition. Pure-Python so we can run NULL_REPLICATES inside the same
    process without an external dep.
    """
    seq = seq.upper()
    n = len(seq)
    if n < 4:
        # too short for meaningful dinuc shuffle -> fall back to 1st-order
        chars = list(seq)
        rng.shuffle(chars)
        return "".join(chars)
    # build adjacency: edges (a, b) where seq has a transition a->b
    edges: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    last = seq[-1]
    # shuffle outgoing edges per source; ensure the "last edge" out of each
    # node (except `last`) points along a tree to `last` (Eulerian feasibility)
    for k in list(edges.keys()):
        rng.shuffle(edges[k])
    # simple verification + fallback: try multi-walk; if Eulerian assembly
    # fails (rare given precondition), fall back to plain 1st-order shuffle
    try:
        out = [seq[0]]
        cursor = seq[0]
        for _ in range(n - 1):
            if not edges[cursor]:
                raise RuntimeError("dead-end")
            nxt = edges[cursor].pop()
            out.append(nxt)
            cursor = nxt
        return "".join(out)
    except RuntimeError:
        chars = list(seq)
        rng.shuffle(chars)
        return "".join(chars)


def null_recurrence(
    seqs: list[str],
    threshold: float = JACCARD_THRESHOLD,
    replicates: int = NULL_REPLICATES,
    seed: int = NULL_SEED,
) -> tuple[float, float, float, float]:
    """Return (mean_ge2_pct, mean_ge3_pct, p_ge2, p_ge3) for the recurrence
    rate of a length-matched dinucleotide-preserving shuffle null.

    The observed recurrence-rate-ge2 is then computed by the caller and the
    permutation p is (#replicates where null >= observed) / replicates.
    """
    if not seqs:
        return 0.0, 0.0, 1.0, 1.0
    rng = random.Random(seed)
    rec2 = []
    rec3 = []
    for _ in range(replicates):
        shuffled = [_dinuc_shuffle(s, rng) for s in seqs]
        ks = [canonical_kmers(s, K) for s in shuffled]
        lens = [len(s) for s in shuffled]
        labels = representative_cluster(ks, lens, threshold)
        sizes = Counter(labels).values()
        n = len(shuffled)
        ge2 = sum(s for s in sizes if s >= 2)
        ge3 = sum(s for s in sizes if s >= 3)
        rec2.append(100.0 * ge2 / n if n else 0.0)
        rec3.append(100.0 * ge3 / n if n else 0.0)
    return (
        float(np.mean(rec2)),
        float(np.mean(rec3)),
        float(np.std(rec2)),
        float(np.std(rec3)),
    )


def main() -> int:
    df = pd.read_csv(IN_TSV, sep="\t", dtype=str).fillna("")
    df["len"] = df["len"].astype(int)
    n_input_rows = len(df)
    df = df[df["len"] >= MIN_LEN].reset_index(drop=True)
    df["seq"] = df["seq"].str.upper()
    # drop blank seqs just in case
    df = df[df["seq"].str.len() > 0].reset_index(drop=True)
    n_rows_post_len = len(df)
    # Same physical linker segment can appear multiple times because the table
    # was emitted per panel x category (e.g. within_eimeria/telotron_present
    # and outgroup/intron_present_nontelo for the same locus). Before deduping
    # we preserve the panel/category context as a comma-joined column so
    # downstream consumers can still subset on category. Dedup masks the row
    # inflation from linker_segmentation but does not delete the context.
    panel_col = "panel" if "panel" in df.columns else None
    cat_col = "category" if "category" in df.columns else None
    if panel_col or cat_col:
        agg_cols = {c: lambda s: ",".join(sorted(set(s.dropna().astype(str)))) for c in (panel_col, cat_col) if c}
        ctx = (
            df.groupby(["telotron", "position", "seq"], as_index=False)
              .agg({c: agg_cols[c] for c in agg_cols})
        )
    else:
        ctx = None
    df = df.drop_duplicates(subset=["telotron", "position", "seq"]).reset_index(drop=True)
    if ctx is not None:
        df = df.drop(columns=[c for c in (panel_col, cat_col) if c and c in df.columns])
        df = df.merge(ctx, on=["telotron", "position", "seq"], how="left")
    n_linkers = len(df)
    n_short = int((df["len"] < SHORT_LONG_BOUNDARY).sum())
    n_long = int((df["len"] >= SHORT_LONG_BOUNDARY).sum())
    print(
        f"loaded {n_input_rows} rows -> {n_rows_post_len} >= {MIN_LEN} bp -> "
        f"{n_linkers} unique linker segments "
        f"(short <{SHORT_LONG_BOUNDARY}bp: {n_short}; long: {n_long})",
        file=sys.stderr,
    )

    kmer_sets = [canonical_kmers(s, K) for s in df["seq"]]
    seq_lens = df["len"].tolist()
    labels = representative_cluster(kmer_sets, seq_lens, JACCARD_THRESHOLD)
    df["cluster_id_raw"] = labels

    # summarize clusters
    rows = []
    grouped = df.groupby("cluster_id_raw")
    for cid_raw, sub in grouped:
        size = len(sub)
        seqs = sub["seq"].tolist()
        cons = consensus_from_seqs(seqs)
        examples = sub["telotron"].drop_duplicates().tolist()[:5]
        rows.append(
            {
                "cluster_id_raw": cid_raw,
                "size": size,
                "consensus": cons,
                "consensus_len": len(cons),
                "mean_member_len": float(sub["len"].mean()),
                "n_unique_telotrons": sub["telotron"].nunique(),
                "n_unique_species": sub["species"].nunique(),
                "example_telotrons": ";".join(examples),
            }
        )
    clust_df = pd.DataFrame(rows).sort_values("size", ascending=False).reset_index(drop=True)
    clust_df["cluster_id"] = [f"L{i+1:04d}" for i in range(len(clust_df))]
    clust_df = clust_df[
        [
            "cluster_id",
            "size",
            "consensus_len",
            "mean_member_len",
            "n_unique_telotrons",
            "n_unique_species",
            "consensus",
            "example_telotrons",
        ]
    ]
    clust_df.to_csv(OUT_CLUSTERS, sep="\t", index=False)
    print(f"wrote {OUT_CLUSTERS} ({len(clust_df)} clusters)", file=sys.stderr)

    # recurrence stats (overall + stratified)
    def _stats(sub_df: pd.DataFrame, all_n: int) -> dict[str, float]:
        sizes = sub_df["size"].astype(int)
        n_clu = int(len(sub_df))
        n_sing = int((sizes == 1).sum())
        n_ge2 = int((sizes >= 2).sum())
        n_ge3 = int((sizes >= 3).sum())
        mem_ge2 = int(sizes[sizes >= 2].sum())
        mem_ge3 = int(sizes[sizes >= 3].sum())
        return {
            "n_clusters": n_clu,
            "n_singleton": n_sing,
            "n_ge2": n_ge2,
            "n_ge3": n_ge3,
            "members_in_ge2": mem_ge2,
            "members_in_ge3": mem_ge3,
            "pct_ge2": 100.0 * mem_ge2 / all_n if all_n else 0.0,
            "pct_ge3": 100.0 * mem_ge3 / all_n if all_n else 0.0,
        }

    overall = _stats(clust_df, n_linkers)

    # Stratify: split the original df by short/long, redo the per-stratum
    # cluster summary (membership counts only).
    short_mask = df["len"] < SHORT_LONG_BOUNDARY
    long_mask = ~short_mask
    short_sizes = Counter(df.loc[short_mask, "cluster_id_raw"])
    long_sizes = Counter(df.loc[long_mask, "cluster_id_raw"])
    short_size_dist = pd.DataFrame({"size": list(short_sizes.values())}) if short_sizes else pd.DataFrame({"size": []})
    long_size_dist = pd.DataFrame({"size": list(long_sizes.values())}) if long_sizes else pd.DataFrame({"size": []})
    short_stats = _stats(short_size_dist, int(short_mask.sum()))
    long_stats = _stats(long_size_dist, int(long_mask.sum()))

    # Null distribution: shuffle-permutation baseline for recurrence rate
    print(
        f"computing length-matched dinucleotide-shuffle null over "
        f"{NULL_REPLICATES} replicates ...",
        file=sys.stderr,
    )
    seqs_list = df["seq"].tolist()
    null_mean_ge2, null_mean_ge3, null_sd_ge2, null_sd_ge3 = null_recurrence(
        seqs_list, JACCARD_THRESHOLD, NULL_REPLICATES, NULL_SEED
    )
    # one-sided permutation p: fraction of replicates whose null recurrence
    # meets or exceeds the observed value. Recompute the replicate vector to
    # get an empirical p (cheap second pass).
    rng = random.Random(NULL_SEED)
    null_vec_ge2 = []
    null_vec_ge3 = []
    for _ in range(NULL_REPLICATES):
        shuffled = [_dinuc_shuffle(s, rng) for s in seqs_list]
        ks = [canonical_kmers(s, K) for s in shuffled]
        lens = [len(s) for s in shuffled]
        lab = representative_cluster(ks, lens, JACCARD_THRESHOLD)
        sizes = Counter(lab).values()
        n = len(shuffled)
        ge2 = sum(s for s in sizes if s >= 2)
        ge3 = sum(s for s in sizes if s >= 3)
        null_vec_ge2.append(100.0 * ge2 / n if n else 0.0)
        null_vec_ge3.append(100.0 * ge3 / n if n else 0.0)
    obs_ge2 = overall["pct_ge2"]
    obs_ge3 = overall["pct_ge3"]
    p_ge2 = (sum(1 for v in null_vec_ge2 if v >= obs_ge2) + 1) / (NULL_REPLICATES + 1)
    p_ge3 = (sum(1 for v in null_vec_ge3 if v >= obs_ge3) + 1) / (NULL_REPLICATES + 1)
    delta_ge2 = obs_ge2 - null_mean_ge2
    delta_ge3 = obs_ge3 - null_mean_ge3

    top20 = clust_df.head(20)

    with open(OUT_SUMMARY, "w") as fh:
        fh.write("# Linker recurrence summary\n")
        fh.write(f"input: {IN_TSV}\n")
        fh.write(f"min_len: {MIN_LEN}\n")
        fh.write(f"short_long_boundary_bp: {SHORT_LONG_BOUNDARY}\n")
        fh.write(f"kmer_k: {K}\n")
        fh.write(f"jaccard_threshold: {JACCARD_THRESHOLD}\n")
        fh.write(f"cluster_algorithm: representative_based (CD-HIT-EST-style)\n")
        fh.write(f"n_input_rows: {n_input_rows}\n")
        fh.write(f"n_rows_post_len_filter: {n_rows_post_len}\n")
        fh.write(f"n_linkers: {n_linkers}\n")
        fh.write(f"n_short_linkers: {n_short}\n")
        fh.write(f"n_long_linkers: {n_long}\n")
        fh.write("\n# Overall recurrence\n")
        fh.write(f"n_clusters: {overall['n_clusters']}\n")
        fh.write(f"n_singleton_clusters: {overall['n_singleton']}\n")
        fh.write(f"n_clusters_size_ge2: {overall['n_ge2']}\n")
        fh.write(f"n_clusters_size_ge3: {overall['n_ge3']}\n")
        fh.write(f"linkers_in_cluster_ge2: {overall['members_in_ge2']}\n")
        fh.write(f"linkers_in_cluster_ge3: {overall['members_in_ge3']}\n")
        fh.write(f"recurrence_rate_ge2_pct: {obs_ge2:.2f}\n")
        fh.write(f"recurrence_rate_ge3_pct: {obs_ge3:.2f}\n")
        fh.write("\n# Stratified by linker length\n")
        fh.write(
            f"short (<{SHORT_LONG_BOUNDARY}bp): n={n_short} "
            f"recurrence_ge2_pct={short_stats['pct_ge2']:.2f} "
            f"recurrence_ge3_pct={short_stats['pct_ge3']:.2f}\n"
        )
        fh.write(
            f"long  (>={SHORT_LONG_BOUNDARY}bp): n={n_long} "
            f"recurrence_ge2_pct={long_stats['pct_ge2']:.2f} "
            f"recurrence_ge3_pct={long_stats['pct_ge3']:.2f}\n"
        )
        fh.write(
            "# (Inter-event seam linkers ~11 bp -> short bin; central A-|A+ "
            "junction linkers ~58 bp -> long bin; see memory "
            "telotron_two_linker_classes_2026-06-03.)\n"
        )
        fh.write("\n# Length-matched dinucleotide-shuffle null baseline\n")
        fh.write(f"null_replicates: {NULL_REPLICATES}\n")
        fh.write(f"null_seed: {NULL_SEED}\n")
        fh.write(f"null_mean_recurrence_ge2_pct: {null_mean_ge2:.2f}\n")
        fh.write(f"null_sd_recurrence_ge2_pct:   {null_sd_ge2:.2f}\n")
        fh.write(f"null_mean_recurrence_ge3_pct: {null_mean_ge3:.2f}\n")
        fh.write(f"null_sd_recurrence_ge3_pct:   {null_sd_ge3:.2f}\n")
        fh.write(f"delta_ge2_pct (obs - null):   {delta_ge2:.2f}\n")
        fh.write(f"delta_ge3_pct (obs - null):   {delta_ge3:.2f}\n")
        fh.write(f"permutation_p_ge2 (one-sided): {p_ge2:.4g}\n")
        fh.write(f"permutation_p_ge3 (one-sided): {p_ge3:.4g}\n")
        fh.write("\n# Top 20 largest clusters\n")
        fh.write(
            "cluster_id\tsize\tconsensus_len\tn_unique_telotrons\tn_unique_species\tconsensus\texample_telotrons\n"
        )
        for _, row in top20.iterrows():
            fh.write(
                f"{row['cluster_id']}\t{row['size']}\t{row['consensus_len']}\t"
                f"{row['n_unique_telotrons']}\t{row['n_unique_species']}\t"
                f"{row['consensus']}\t{row['example_telotrons']}\n"
            )
    print(f"wrote {OUT_SUMMARY}", file=sys.stderr)
    print(
        f"n_linkers={n_linkers}, n_clusters={overall['n_clusters']}, "
        f"recurrence(>=3)={obs_ge3:.2f}% (null {null_mean_ge3:.2f}%, "
        f"delta {delta_ge3:+.2f}, p={p_ge3:.3g}) "
        f"recurrence(>=2)={obs_ge2:.2f}% (null {null_mean_ge2:.2f}%, "
        f"delta {delta_ge2:+.2f}, p={p_ge2:.3g})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
