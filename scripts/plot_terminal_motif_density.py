#!/usr/bin/env python3
"""Per-species figures: rolling-average density of each genome's identified terminal
telomere motif across its top-N longest contigs. Telotron loci overlaid as light marks."""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(it=None, **kw):
        return it if it is not None else iter(())
    tqdm.write = print  # type: ignore[attr-defined]

from _common import rotations, all_variants, find_genome_fasta
from telomere_mask import mask_telomere_fragments


def infer_terminal_motif(fasta_path, candidate_motifs, end_window=1000):
    """Scan first+last end_window bp of every contig for 5-12mer frequencies.
    For each candidate, score = (hits at k=len(candidate)) / (total windows at that k).
    This normalises across different repeat unit lengths so a 5-mer doesn't beat a
    6-mer just because there are more 5-mer windows.
    Returns (best_motif, kmer_report_str).
    """
    import gzip
    from collections import Counter, defaultdict

    # Build lookup: kmer_variant -> candidate motif (variant length == candidate length)
    lookup = {}
    for motif in candidate_motifs:
        for v in all_variants(motif):
            if v not in lookup:
                lookup[v] = motif

    # Stream FASTA; only accumulate terminal end_window characters per contig.
    # Use a fixed-size deque to avoid buffering multi-MB chromosomes.
    from collections import deque

    opener = gzip.open if fasta_path.endswith(".gz") else open
    terminal_seqs = []

    def emit(head_buf, tail_buf, length):
        w = min(end_window, length // 2)
        if w <= 0:
            return
        head = "".join(head_buf)[:w].upper()
        tail = "".join(tail_buf)[-w:].upper()
        if head:
            terminal_seqs.append(head)
        if tail and tail != head:
            terminal_seqs.append(tail)

    with opener(fasta_path, "rt") as fh:
        head_buf, tail_buf, length = [], deque(), 0
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if length:
                    emit(head_buf, tail_buf, length)
                head_buf, tail_buf, length = [], deque(), 0
            else:
                length += len(line)
                # collect head until we have end_window chars
                if sum(len(s) for s in head_buf) < end_window:
                    head_buf.append(line)
                # always update tail (sliding buffer of last end_window chars)
                tail_buf.append(line)
                while sum(len(s) for s in tail_buf) > end_window + 200:
                    tail_buf.popleft()
        if length:
            emit(head_buf, tail_buf, length)

    if not terminal_seqs:
        return "", "no terminal sequences"

    # Count all 5-12mers; track raw counts for reporting and per-k totals for normalisation.
    k_min, k_max = 5, 13
    raw_counts = Counter()          # (k, kmer) -> count
    cand_hits = Counter()           # candidate -> hits at k=len(candidate)
    total_at_k = Counter()          # k -> total valid windows

    for seq in terminal_seqs:
        n = len(seq)
        for k in range(k_min, k_max):
            for i in range(n - k + 1):
                kmer = seq[i:i + k]
                if all(c in "ACGT" for c in kmer):
                    raw_counts[(k, kmer)] += 1
                    total_at_k[k] += 1
                    cand = lookup.get(kmer)
                    if cand:
                        cand_hits[cand] += 1

    # ---- 1st-order Markov background from telomere-MASKED terminal sequences ----
    # Uniform 4^k background overweights long T-rich motifs on AT-rich (e.g. 80% AT
    # Eimeria) genomes; per MEMORY eimeria_subtelomere_rlfs_2026-05-15, any AT-context
    # composition test on telomere-bearing sequence must mask telomere rotations + revcomp
    # + 1-mm variants first. Build mono- and di-nucleotide frequencies from the masked
    # terminal windows; expected fraction for a k-mer is P(b0) * prod P(b_i | b_{i-1}).
    mono = Counter()
    di = Counter()
    for seq in terminal_seqs:
        masked = mask_telomere_fragments(seq)
        prev = None
        for c in masked:
            if c not in "ACGT":
                prev = None
                continue
            mono[c] += 1
            if prev is not None:
                di[(prev, c)] += 1
            prev = c
    mono_total = sum(mono.values())
    if mono_total == 0:
        # Fallback to uniform if masking erased everything (pure-telomere terminus).
        p_mono = {c: 0.25 for c in "ACGT"}
        p_cond = {(a, b): 0.25 for a in "ACGT" for b in "ACGT"}
    else:
        # Laplace pseudocount so absent bases / dinucleotides don't zero-out P.
        p_mono = {c: (mono.get(c, 0) + 1) / (mono_total + 4) for c in "ACGT"}
        p_cond = {}
        for a in "ACGT":
            row_total = sum(di.get((a, b), 0) for b in "ACGT")
            for b in "ACGT":
                p_cond[(a, b)] = (di.get((a, b), 0) + 1) / (row_total + 4)

    def _markov_prob(kmer):
        p = p_mono[kmer[0]]
        for i in range(1, len(kmer)):
            p *= p_cond[(kmer[i - 1], kmer[i])]
        return p

    # Normalise: fold-enrichment = observed fraction / expected Markov fraction
    # (summed over all rotation + revcomp variants of the candidate motif).
    # Score each candidate by: enrichment × hits_per_window.
    # - enrichment = observed_fraction / expected_markov_fraction (k-length + AT-context aware)
    # - hits_per_window = hits / n_terminal_windows (rewards consistent signal across contigs
    #   vs. isolated spurious matches of rare long k-mers)
    # Both filters must pass: enrichment ≥ MIN_ENRICHMENT and hits ≥ MIN_HITS.
    MIN_ENRICHMENT = 3.0
    MIN_HITS = 20
    n_windows = max(1, len(terminal_seqs))
    cand_enrichment = {}
    cand_score = {}
    for cand, hits in cand_hits.items():
        if hits < MIN_HITS:
            continue
        k = len(cand)
        expected_frac = sum(_markov_prob(v) for v in all_variants(cand))
        observed_frac = hits / max(1, total_at_k[k])
        enrich = observed_frac / expected_frac if expected_frac else 0.0
        if enrich < MIN_ENRICHMENT:
            continue
        cand_enrichment[cand] = enrich
        cand_score[cand] = enrich * (hits / n_windows)

    best = max(cand_score, key=cand_score.get) if cand_score else ""

    # Report top 3 k-mers per length
    lines = []
    for k in range(k_min, k_max):
        top = sorted(((cnt, km) for (klen, km), cnt in raw_counts.items() if klen == k),
                     reverse=True)[:3]
        if top:
            lines.append(f"  k={k}: " + "  ".join(f"{km}={cnt}" for cnt, km in top))
    if cand_enrichment:
        lines.append("  candidate enrichment (enrichment × n_hits): " +
                     "  ".join(f"{c}={cand_enrichment[c]:.1f}x(n={cand_hits[c]})"
                               for c in sorted(cand_enrichment, key=lambda x: -cand_score.get(x, 0))[:5]))
    report = "\n".join(lines) if lines else "  (none)"

    return best, report


def contig_lengths(fasta_path):
    r = subprocess.run(["seqkit", "fx2tab", "-j", "1", "-i", "-l", "-n", fasta_path],
                       capture_output=True, text=True, check=True)
    out = {}
    for line in r.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        out[parts[0]] = int(parts[-1])
    return out


def locate_motif(fasta_path, motif, bed_path):
    """seqkit hits for all rotations of `motif` on both strands."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as patt:
        for rot in rotations(motif):
            patt.write(f">{motif}\n{rot}\n")
        patt_path = patt.name
    try:
        with open(bed_path, "w") as out:
            subprocess.run(["seqkit", "locate", "-j", "1", "--bed", "--ignore-case",
                            "--pattern-file", patt_path, fasta_path],
                           check=True, stdout=out, stderr=subprocess.DEVNULL)
    finally:
        os.unlink(patt_path)


def per_contig_params(length):
    """Adaptive bin size + smoothing window so a ~500 bp telomere stays visible.
    Rule of thumb: target ~bin = telomere_width / 5, smoothing ~ bin × 3-10."""
    if length < 50_000:
        return 20, 3       # 60 bp smoothing
    if length < 500_000:
        return 100, 5      # 500 bp smoothing
    if length < 5_000_000:
        return 200, 5      # 1 kb smoothing
    return 500, 10         # 5 kb smoothing


def plot_one(args):
    (gid, organism, motif, motif_src, fasta, loci_df, outpath,
     top_n, min_contig_len) = args
    if motif:
        with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as bed:
            bed_path = bed.name
        try:
            locate_motif(fasta, motif, bed_path)
            if os.path.getsize(bed_path):
                hits = pd.read_csv(bed_path, sep="\t", header=None,
                                   names=["seqid", "start", "end", "name", "score", "strand"])
            else:
                hits = pd.DataFrame(columns=["seqid", "start", "end", "name", "score", "strand"])
        finally:
            os.unlink(bed_path)
    else:
        hits = pd.DataFrame(columns=["seqid", "start", "end", "name", "score", "strand"])

    lens = contig_lengths(fasta)
    contigs = [(s, L) for s, L in lens.items() if L >= min_contig_len]
    contigs.sort(key=lambda kv: -kv[1])
    contigs = contigs[:top_n]
    if not contigs:
        return f"skip {gid}: no contigs ≥ {min_contig_len} bp"

    n = len(contigs)

    # Pre-group hits by seqid once (avoid per-contig pandas filter)
    hits_by_seq = {s: g for s, g in hits.groupby("seqid")} if not hits.empty else {}
    loci_by_seq = {s: g for s, g in loci_df.groupby("seqid")} if loci_df is not None and not loci_df.empty else {}

    DISPLAY_BINS = 1200            # ~panel-pixel-width target
    SMOOTH_KB = 2                  # rolling-mean physical window (kb), constant across panels

    def maxpool(arr, n_out):
        """Max-pool 1-D array down to n_out bins. Preserves single-bin peaks."""
        n = len(arr)
        if n <= n_out:
            return arr
        factor = (n + n_out - 1) // n_out
        pad = factor * n_out - n
        if pad > 0:
            arr = np.concatenate([arr, np.zeros(pad, dtype=arr.dtype)])
        return arr.reshape(n_out, factor).max(axis=1)

    # ---- Pass 1: compute density (hits/kb) for every contig, find global max ----
    panels = []
    for seqid, length in contigs:
        bs, _ = per_contig_params(length)
        n_bins = max(1, length // bs + 1)
        counts = np.zeros(n_bins, dtype=float)
        h = hits_by_seq.get(seqid)
        if h is not None and not h.empty:
            idx = (h.start.values // bs).astype(int)
            idx = idx[(idx >= 0) & (idx < n_bins)]
            np.add.at(counts, idx, 1)

        # Convert to density (hits/kb) so panels with different bin sizes are comparable.
        density = counts / (bs / 1000.0)
        win_bins = max(1, int(round(SMOOTH_KB * 1000.0 / bs)))
        smooth = pd.Series(density).rolling(window=win_bins, center=True, min_periods=1).mean().values
        display = maxpool(smooth, DISPLAY_BINS)
        panels.append((seqid, length, display))

    global_ymax = max(max(d.max(), 1.0) for _, _, d in panels) * 1.1
    max_length = max(length for _, length, _ in panels)

    # ---- Pass 2: genome-browser layout — physical panel width ∝ contig length ----
    # Each panel's physical width in figure-fraction = AVAIL_W × (length / max_length),
    # so 1 inch represents the same number of bp in every panel.
    FIG_W = 10.0          # inches
    PANEL_H = 0.55        # inches per panel
    fig_h = max(2.0, n * PANEL_H + 0.5)
    fig = plt.figure(figsize=(FIG_W, fig_h))

    LEFT_IN = 2.1         # inches reserved for ylabel
    RIGHT_IN = 0.15       # right margin
    AVAIL_IN = FIG_W - LEFT_IN - RIGHT_IN   # inches available for data at max_length
    LEFT = LEFT_IN / FIG_W
    AVAIL_W = AVAIL_IN / FIG_W
    DATA_H = PANEL_H * 0.72 / fig_h        # panel data area height in fig fraction
    GAP_H = PANEL_H * 0.28 / fig_h         # gap between panels

    axes = []
    for i, (seqid, length, display) in enumerate(panels):
        w = AVAIL_W * (length / max_length)
        b = 1.0 - (i + 1) * (PANEL_H / fig_h) + GAP_H / 2
        ax = fig.add_axes([LEFT, b, w, DATA_H])
        axes.append(ax)

    for ax, (seqid, length, display) in zip(axes, panels):
        x_mb = np.linspace(0, length / 1e6, len(display))
        ax.fill_between(x_mb, display, color="#A84B2F", alpha=0.85, linewidth=0)

        l = loci_by_seq.get(seqid)
        if l is not None and not l.empty:
            ax.scatter(l.start.values / 1e6, [global_ymax * 0.95] * len(l),
                       marker="v", s=8, color="#003366", alpha=0.85, linewidths=0,
                       clip_on=False)

        if length >= 1_000_000:
            len_str = f"{length / 1e6:.1f} Mb"
        elif length >= 1_000:
            len_str = f"{length / 1e3:.0f} kb"
        else:
            len_str = f"{length} bp"
        ax.set_ylabel(f"{seqid[:18]}\n{len_str}", rotation=0, ha="right", va="center", fontsize=7)
        ax.set_xlim(0, length / 1e6)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.set_yticks([])
        ax.tick_params(which="both", left=False)
        # set_ylim last so nothing overrides it; disable autoscale explicitly
        ax.set_autoscaley_on(False)
        ax.set_ylim(0, global_ymax)

    axes[-1].set_xlabel("Genomic position (Mb)")
    if not motif:
        motif_label = "none"
    elif motif_src == "kmer":
        motif_label = f"{motif} (k-mer inferred)"
    else:
        motif_label = motif
    fig.suptitle(f"{organism}  ({gid})  —  motif: {motif_label}  "
                 f"—  rolling-mean hits/kb (shared y-axis, {SMOOTH_KB} kb window)  "
                 f"—  ▼ telotron candidates",
                 fontsize=10)
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return f"wrote {outpath}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-n", type=int, default=20,
                    help="Plot the N longest contigs per species (one panel each).")
    ap.add_argument("--min-contig-length", type=int, default=10000)
    ap.add_argument("--motifs", default="",
                    help="Comma-separated candidate telomere motifs (G-rich strand). "
                         "Used to infer terminal motif when contig-end scan found none.")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    candidate_motifs = [m.strip() for m in args.motifs.split(",") if m.strip()] if args.motifs else []

    summary = pd.read_csv(args.summary, sep="\t")
    try:
        loci = pd.read_csv(args.loci, sep="\t")
    except (pd.errors.EmptyDataError, FileNotFoundError):
        loci = pd.DataFrame()

    os.makedirs(args.outdir, exist_ok=True)
    jobs = []
    for r in summary.itertuples(index=False):
        if r.status != "OK":
            continue
        motif = "" if pd.isna(r.terminal_motif) else r.terminal_motif
        motif_src = "terminal" if motif else ""
        try:
            fa = find_genome_fasta(r.genome_id, args.refseq_dir, args.tara_dir)
        except FileNotFoundError as e:
            sys.stderr.write(f"skip {r.genome_id}: {e}\n")
            continue
        if not motif and candidate_motifs:
            inferred, report = infer_terminal_motif(fa, candidate_motifs)
            if inferred:
                motif = inferred
                motif_src = "kmer"
                sys.stderr.write(f"{r.genome_id}: inferred motif={motif} from k-mer scan\n{report}\n")
            else:
                sys.stderr.write(f"{r.genome_id}: no motif inferred\n{report}\n")
        spc_loci = loci[loci.genome_id == r.genome_id] if len(loci) else None
        safe = re.sub(r"[^A-Za-z0-9_]+", "_", f"{r.genome_id}_{r.organism}")
        out = os.path.join(args.outdir, f"{safe}.png")
        jobs.append((r.genome_id, r.organism, motif, motif_src, fa, spc_loci, out,
                     args.top_n, args.min_contig_length))

    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        futures = [ex.submit(plot_one, j) for j in jobs]
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc="plotting species", unit="species"):
            tqdm.write(fut.result())


if __name__ == "__main__":
    main()
