#!/usr/bin/env python3
"""Boundary k-mer enrichment, distance-to-contig-end, and architecture summary.

**Rotation-trap note:** `boundary_kmers()` counts k-mers at the donor/acceptor
edges of telotrons vs non-telotron introns. Telotron edges are BY CONSTRUCTION
rotations of the defining telomere motif, so enrichment is definitional. The
output is flagged `circular_by_construction=True` and consumers should treat
the Fisher p-values as descriptive, not inferential. Telomere-fragment masking
(`telomere_mask.mask_telomere_fragments`) does not apply here because the
k-mer counts themselves ARE the signal — masking would delete what we're
counting. The non-circular question (register bias beyond splice-signal
selection) is a within-array null handled by the splice-site-gate analysis."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from _common import tqdm

try:
    from scipy.stats import fisher_exact, binomtest
except Exception:
    fisher_exact = None
    binomtest = None


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
    """Per (genome, motif), top 10 k-mers at the donor/acceptor edges vs non-telotron introns.

    Fisher's exact p + BH q are reported for continuity with prior versions,
    but the comparison is DEFINITIONALLY CIRCULAR: telotron boundary k-mers
    are rotations of the telomere repeat that defines a telotron, so
    enrichment over non-telomeric introns is expected by construction. Rows
    with fold_enrichment='inf' come from k-mers absent from every control
    intron (they are all definitional rotations of the motif); every table row
    is flagged with `circular_by_construction=True` and downstream consumers
    should treat the Fisher p as descriptive, not inferential. The
    non-circular question (register bias beyond splice-signal selection) is a
    within-array null handled by the splice-site-gate analysis."""
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

    df = pd.DataFrame(rows, columns=[
        "genome_id", "motif", "boundary", "kmer",
        "telotron_count", "telotron_freq",
        "control_count", "control_freq", "fold_enrichment",
        "fisher_p", "bh_q",
    ])
    # Flag every row: this table is definitionally circular (see docstring).
    # Downstream consumers should treat p/q values as descriptive, not inferential.
    df["circular_by_construction"] = True
    df.to_csv(out_path, sep="\t", index=False)


def _distance_one(args):
    """Per-genome subtelomeric-distance test against a WITHIN-CONTIG null.

    The previous pooled Mann-Whitney treated co-contig loci as independent
    (pseudoreplication) and compared against a duplicated, sampled control set.
    Instead, each telotron is judged against the distance-to-end expected if an
    intron of the same length were placed uniformly at random on ITS OWN contig:
    distance = min(U, free-U), U ~ Unif(0, free), free = contig_len - intron_len.
    This removes cross-locus dependence and control duplication, and conditions
    on contig length (short contigs mechanically force small distances).
    """
    genome_id, group, capped_only, n_perm, seed = args
    if capped_only and "contig_both_ends_capped" in group.columns:
        keep = group["contig_both_ends_capped"].astype(str).str.lower().isin(["true", "1"])
        group = group[keep]
    n = len(group)
    if not n:
        return [genome_id, 0, 0, "", "", "", "", ""]
    obs = group.distance_to_end.astype(float).to_numpy()
    free = (group.contig_len.astype(float)
            - group.intron_len.astype(float)).clip(lower=1.0).to_numpy()
    obs_med = float(np.median(obs))
    # Sign test: is each telotron closer to its contig end than that contig's
    # null median (free/4)?  Binomial vs 0.5, one-sided (closer than chance).
    closer = int(np.count_nonzero(obs < free / 4.0))
    sign_p = (binomtest(closer, n, 0.5, alternative="greater").pvalue
              if binomtest else "")
    # Monte-Carlo: null distribution of the median distance under random placement.
    rng = np.random.default_rng(seed)
    null_meds = np.empty(n_perm)
    for i in range(n_perm):
        u = rng.uniform(0.0, free)
        null_meds[i] = np.median(np.minimum(u, free - u))
    mc_p = (np.count_nonzero(null_meds <= obs_med) + 1) / (n_perm + 1)
    return [genome_id, n, int(group.seqid.nunique()), obs_med,
            float(np.median(null_meds)), closer / n, sign_p, mc_p]


def distance_to_end(final, out_path, threads=1, capped_only=True, n_perm=2000, seed=1):
    """Per-genome subtelomeric-clustering test vs a within-contig random-placement
    null. Emits BOTH the capped-only pass (telomere-capped contigs, defensible
    for chromosome-level assemblies) and an `_all_contigs` companion (every
    contig, defensible for MAGs where "chromosome end" isn't meaningful — the
    capped-only test drops most Tara MAG loci and gives a high-variance zero)."""
    def _run(_capped, suffix):
        tasks = [(gid, group, _capped, n_perm, seed)
                 for gid, group in final.groupby("genome_id")]
        if threads > 1 and len(tasks) > 1:
            with ProcessPoolExecutor(max_workers=threads) as ex:
                rows = list(tqdm(ex.map(_distance_one, tasks), total=len(tasks),
                                 desc=f"distance-to-end{suffix}", unit="genome"))
        else:
            rows = [_distance_one(t) for t in tqdm(tasks, desc=f"distance-to-end{suffix}", unit="genome")]
        df = pd.DataFrame(rows, columns=[
            "genome_id", "telotron_n", "n_capped_contigs",
            "median_telotron_distance_to_end", "median_null_distance_to_end",
            "frac_closer_than_contig_null", "sign_test_p", "mc_permutation_p",
        ])
        df["mode"] = "capped_only" if _capped else "all_contigs"
        return df
    pd.concat([_run(True, " (capped)"), _run(False, " (all)")], ignore_index=True) \
        .to_csv(out_path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True)
    ap.add_argument("--introns", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--boundary-kmers", required=True)
    ap.add_argument("--distance", required=True)
    ap.add_argument("--architecture", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--all-contigs", action="store_true",
                    help="distance test: use all contigs, not just fully "
                         "telomere-capped ones (default: capped contigs only).")
    args = ap.parse_args()

    final = pd.read_csv(args.final, sep="\t")
    introns = pd.read_csv(args.introns, sep="\t", low_memory=False)

    boundary_kmers(final, introns, args.boundary_kmers, threads=args.threads)
    distance_to_end(final, args.distance, threads=args.threads,
                    capped_only=not args.all_contigs)
    (final.groupby(["genome_id", "orientation", "splice_class"])
          .size().reset_index(name="n")
          .to_csv(args.architecture, sep="\t", index=False))


if __name__ == "__main__":
    main()
