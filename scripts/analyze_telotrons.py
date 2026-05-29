#!/usr/bin/env python3
"""Boundary k-mer enrichment, distance-to-contig-end, and architecture summary."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from _common import tqdm

try:
    from scipy.stats import mannwhitneyu, fisher_exact
except Exception:
    mannwhitneyu = None
    fisher_exact = None


LOCUS_KEY = ["genome_id", "seqid", "start", "end", "strand"]


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values for a list of p-values (None-safe)."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    m = len(idx)
    q = [None] * len(pvals)
    if not m:
        return q
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        q[i] = prev = val
    return q


def split_controls(final, introns):
    """Introns minus the final telotron loci."""
    keys = set(map(tuple, final[LOCUS_KEY].values))
    return introns[~introns[LOCUS_KEY].apply(tuple, axis=1).isin(keys)]


def _boundary_one(args):
    """Per-(genome, motif) boundary k-mer rows. Runs in a worker process."""
    genome_id, motif, group, ctrl = args
    k = len(str(motif))
    positions = [
        ("donor_inside",    "first40", lambda s: s[2:2 + k]),
        ("acceptor_inside", "last40",  lambda s: s[-(k + 2):-2]),
    ]
    rows = []
    for name, column, slicer in positions:
        telo_counts = Counter(group[column].map(str).map(slicer))
        ctrl_counts = Counter(ctrl[column].map(str).map(slicer))
        telo_total = sum(telo_counts.values()) or 1
        ctrl_total = sum(ctrl_counts.values())
        for mer, n in telo_counts.most_common(10):
            cf = ctrl_counts.get(mer, 0)
            tf = n / telo_total
            cfq = cf / ctrl_total if ctrl_total else 0
            fold = (tf / cfq) if cfq else "inf"
            # 2x2 Fisher's exact: this k-mer vs all others, telotron vs control.
            p = ""
            if fisher_exact and ctrl_total:
                p = fisher_exact(
                    [[n, telo_total - n], [cf, ctrl_total - cf]],
                    alternative="greater",
                ).pvalue
            rows.append([genome_id, motif, name, mer, n, tf, cf, cfq, fold, p])
    return rows


def boundary_kmers(final, introns, out_path, threads=1):
    """Per (genome, motif), top 10 k-mers at the donor/acceptor edges vs non-telotron introns."""
    controls = split_controls(final, introns)
    # Pre-slice control set per genome so each worker only gets the data it needs.
    ctrl_by_genome = {gid: g for gid, g in controls.groupby("genome_id")}
    tasks = [
        (gid, motif, group, ctrl_by_genome.get(gid, controls.iloc[0:0]))
        for (gid, motif), group in final.groupby(["genome_id", "motif"])
    ]
    if threads > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            results = list(tqdm(ex.map(_boundary_one, tasks), total=len(tasks),
                                desc="boundary k-mers", unit="grp"))
    else:
        results = [_boundary_one(t) for t in tqdm(tasks, desc="boundary k-mers", unit="grp")]
    rows = [r for sub in results for r in sub]

    # Benjamini-Hochberg FDR across every k-mer tested (the raw fold-enrichment
    # alone has no significance; thousands of k-mers are compared).
    qvals = bh_fdr([r[9] if r[9] != "" else None for r in rows])
    for r, q in zip(rows, qvals):
        r.append(q if q is not None else "")

    pd.DataFrame(rows, columns=[
        "genome_id", "motif", "boundary", "kmer",
        "telotron_count", "telotron_freq",
        "control_count", "control_freq", "fold_enrichment",
        "fisher_p", "bh_q",
    ]).to_csv(out_path, sep="\t", index=False)


def _distance_one(args):
    """Per-genome distance-to-end stats. Runs in a worker process."""
    genome_id, group, pool = args
    picks = []
    for r in group.itertuples(index=False):
        lo = max(1, 0.8 * float(r.intron_len))
        hi = 1.25 * float(r.intron_len)
        same_contig = pool[(pool.seqid == r.seqid)
                           & (pool.intron_len >= lo) & (pool.intron_len <= hi)]
        candidates = same_contig if len(same_contig) >= 5 else \
            pool[(pool.intron_len >= lo) & (pool.intron_len <= hi)]
        if len(candidates):
            picks.append(candidates.sample(min(10, len(candidates)),
                                           random_state=int(r.start) % 1000003))
    # De-duplicate: the same control intron can be length-matched to several
    # telotrons; counting it once keeps the Mann-Whitney sample independent.
    ctrl = (pd.concat(picks).drop_duplicates(subset=["seqid", "start", "end"])
            if picks else pool.head(0))
    telo_dist = group.distance_to_end.astype(float)
    ctrl_dist = ctrl.distance_to_end.astype(float) if len(ctrl) else pd.Series(dtype=float)
    pval = (mannwhitneyu(telo_dist, ctrl_dist, alternative="two-sided").pvalue
            if mannwhitneyu and len(telo_dist) and len(ctrl_dist) else "")
    return [
        genome_id, len(group), len(ctrl_dist),
        telo_dist.median() if len(telo_dist) else "",
        ctrl_dist.median() if len(ctrl_dist) else "",
        pval,
    ]


def distance_to_end(final, introns, out_path, threads=1):
    """Compare telotron distance-to-contig-end against length-matched non-telotron introns."""
    controls = split_controls(final, introns)
    ctrl_by_genome = {gid: g for gid, g in controls.groupby("genome_id")}
    tasks = [
        (gid, group, ctrl_by_genome.get(gid, controls.iloc[0:0]))
        for gid, group in final.groupby("genome_id")
    ]
    if threads > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            rows = list(tqdm(ex.map(_distance_one, tasks), total=len(tasks),
                             desc="distance-to-end", unit="genome"))
    else:
        rows = [_distance_one(t) for t in tqdm(tasks, desc="distance-to-end", unit="genome")]

    pd.DataFrame(rows, columns=[
        "genome_id", "telotron_n", "length_matched_control_n",
        "median_telotron_distance_to_end", "median_control_distance_to_end",
        "mannwhitney_p",
    ]).to_csv(out_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True)
    ap.add_argument("--introns", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--boundary-kmers", required=True)
    ap.add_argument("--distance", required=True)
    ap.add_argument("--architecture", required=True)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    final = pd.read_csv(args.final, sep="\t")
    introns = pd.read_csv(args.introns, sep="\t", low_memory=False)

    boundary_kmers(final, introns, args.boundary_kmers, threads=args.threads)
    distance_to_end(final, introns, args.distance, threads=args.threads)
    (final.groupby(["genome_id", "orientation", "splice_class"])
          .size().reset_index(name="n")
          .to_csv(args.architecture, sep="\t", index=False))


if __name__ == "__main__":
    main()
