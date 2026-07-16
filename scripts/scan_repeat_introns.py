#!/usr/bin/env python3
"""Broad tandem-repeat scan of introns.

Where scan_telotrons.py restricts the search to the configured telomere
motifs, this script finds introns dominated by ANY short tandem repeat unit
(k = 4-10 by default) and reports:
  - dominant_unit         the observed k-mer with the largest tandem coverage
  - dominant_canonical    lex-min rotation of unit (and of its revcomp) so
                          all rotations/RCs collapse to one identifier
  - unit_len              k
  - unit_frac             matches * k / intron_len  ~ fraction covered
  - telomere_match        True if canonical matches a config telomere motif
  - telomere_match_name   the matching config motif (e.g. TTAGGG), else ""

Output is intended for MANUAL CURATION — no filtering by "known telomere"
happens here. Use --min-frac to gate the noise floor (default 0.5).

Wildcarded per-genome mode: --single-genome <gid> produces just that
genome's rows, so the Snakefile can fan out one sbatch job per genome.
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import rc, rotations, read_fasta, find_genome_files, open_maybe_gz

# scan_telotrons contains the gt-driven intron extractor; reuse it verbatim
# so both scans agree on which introns exist for a given GFF (silent
# gt-truncation trap already handled there).
from scan_telotrons import run_gt_introns, load_introns


K_MIN_DEFAULT = 4
K_MAX_DEFAULT = 10
MIN_FRAC_DEFAULT = 0.5
MIN_INTRON_LEN_DEFAULT = 30


def canonical_unit(unit: str) -> str:
    """Lex-min rotation of unit and of its reverse complement.

    Collapses every rotation and both strands of the same tandem unit to
    a single identifier — TTAGGG, TAGGGT, GGGTTA, CCCTAA all become the
    same string.
    """
    u = unit.upper()
    both = rotations(u) + rotations(rc(u))
    return min(both)


def dominant_repeat(seq: str, k_min: int, k_max: int) -> tuple[str, str, int, float]:
    """Return (unit, canonical_unit, k, frac) for the k-mer with the
    largest tandem coverage in `seq`.

    coverage = (# distinct sliding-window positions matching the dominant
    k-mer) * k / len(seq). For a perfect (TTAGGG)n array this returns ~1.0;
    for a random k-mer occurring twice in a 60 bp intron this returns 0.2.
    """
    n = len(seq)
    if n < 2 * k_min:
        return "", "", 0, 0.0
    best = (0.0, "", "", 0)
    for k in range(k_min, k_max + 1):
        if n < 2 * k:
            continue
        c: Counter = Counter()
        # Sliding-window counts; every position contributes once.
        for i in range(n - k + 1):
            c[seq[i:i + k]] += 1
        # Reject homopolymer/dinucleotide runs — 'A'*k, 'AT'*k/2, etc. hit
        # near-full coverage on any AT-rich intron and would drown out
        # real tandem repeats. Flag as "low_complexity" in the output.
        unit, cnt = c.most_common(1)[0]
        if len(set(unit)) == 1:
            continue                # homopolymer
        frac = min(cnt * k, n) / n
        if frac > best[0]:
            best = (frac, unit, canonical_unit(unit), k)
    return best[1], best[2], best[3], best[0]


def load_telomere_lookup(motifs_csv: str) -> dict[str, str]:
    """Map canonical(motif) -> input motif for every configured telomere motif."""
    out: dict[str, str] = {}
    for m in (motifs_csv or "").split(","):
        m = m.strip().upper()
        if not m:
            continue
        out.setdefault(canonical_unit(m), m)
    return out


def revcomp_if_neg(seq: str, strand: str) -> str:
    return rc(seq) if strand == "-" else seq


def scan_one(args):
    row, opts = args
    gid, org, group, source = row
    try:
        fa, gff = find_genome_files(gid, opts["refseq_dir"], opts["tara_dir"])
        seqs = read_fasta(fa)
    except FileNotFoundError:
        return {"error": f"missing fasta or gff for {gid}", "rows": []}

    with tempfile.TemporaryDirectory(prefix=f"repeat_{gid}_") as tmp:
        gff_with_introns = os.path.join(tmp, "with_introns.gff3")
        try:
            run_gt_introns(gff, gff_with_introns)
        except Exception as e:
            return {"error": f"gt failed for {gid}: {e}", "rows": []}
        introns = load_introns(gff_with_introns)

    if introns.empty:
        return {"error": "", "rows": []}

    telomere_lookup = opts["telomere_lookup"]
    k_min = opts["k_min"]
    k_max = opts["k_max"]
    min_frac = opts["min_frac"]
    min_intron_len = opts["min_intron_len"]

    rows = []
    for r in introns.itertuples(index=False):
        # intron_len already computed in load_introns; recompute defensively.
        s0 = int(r.start) - 1                              # GFF 1-based inclusive → python
        e0 = int(r.end)
        contig = seqs.get(r.seqid)
        if contig is None or e0 > len(contig):
            continue
        ilen = e0 - s0
        if ilen < min_intron_len:
            continue
        raw = contig[s0:e0]
        # Repeats are reported on the transcript (mRNA) strand so a downstream
        # canonical-unit match is orientation-consistent with splice signals.
        seq = revcomp_if_neg(raw, r.strand)
        unit, canon, k, frac = dominant_repeat(seq, k_min, k_max)
        if not unit or frac < min_frac:
            continue
        match_name = telomere_lookup.get(canon, "")
        rows.append([
            gid, org, group, source,
            r.seqid, s0 + 1, e0, r.strand, r.tx_id, r.gene_id, r.intron_index, ilen,
            unit, canon, k, round(frac, 4),
            bool(match_name), match_name,
        ])
    return {"error": "", "rows": rows}


COLUMNS = [
    "genome_id", "organism", "group", "source",
    "seqid", "start", "end", "strand", "tx_id", "gene_id", "intron_index", "intron_len",
    "dominant_unit", "dominant_canonical", "unit_len", "unit_frac",
    "telomere_match", "telomere_match_name",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--single-genome", default="",
                    help="Restrict to this genome_id (slurm fanout).")
    ap.add_argument("--telomere-motifs", default="",
                    help="comma-separated config telomere motifs — matches "
                         "flag telomere_match=True in the output (does NOT "
                         "restrict the search)")
    ap.add_argument("--k-min", type=int, default=K_MIN_DEFAULT)
    ap.add_argument("--k-max", type=int, default=K_MAX_DEFAULT)
    ap.add_argument("--min-frac", type=float, default=MIN_FRAC_DEFAULT,
                    help="only emit introns whose dominant unit covers this "
                         "fraction of the intron (default 0.5)")
    ap.add_argument("--min-intron-len", type=int, default=MIN_INTRON_LEN_DEFAULT)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    rows = []
    with open_maybe_gz(args.manifest) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if args.single_genome and r["genome_id"] != args.single_genome:
                continue
            rows.append([r["genome_id"], r["organism"], r["group"], r["source"]])
    if args.single_genome and not rows:
        sys.exit(f"--single-genome={args.single_genome!r} not present in manifest {args.manifest}")
    if not rows:
        sys.exit(f"no rows in manifest {args.manifest}")

    opts = {
        "refseq_dir": args.refseq_dir,
        "tara_dir": args.tara_dir,
        "telomere_lookup": load_telomere_lookup(args.telomere_motifs),
        "k_min": args.k_min,
        "k_max": args.k_max,
        "min_frac": args.min_frac,
        "min_intron_len": args.min_intron_len,
    }

    tasks = [(row, opts) for row in rows]
    if args.threads > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            results = list(ex.map(scan_one, tasks))
    else:
        results = [scan_one(t) for t in tasks]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n_err = 0
    n_rows = 0
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(COLUMNS)
        for res in results:
            if res["error"]:
                n_err += 1
                print(f"[repeat_scan] {res['error']}", file=sys.stderr)
            for row in res["rows"]:
                w.writerow(row)
                n_rows += 1
    print(f"[repeat_scan] wrote {n_rows} rows across {len(rows)} genomes "
          f"({n_err} errors) → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
