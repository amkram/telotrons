#!/usr/bin/env python3
"""Broad tandem-repeat scan of introns via ULTRA.

For every annotated intron this script runs ULTRA (Olson & Wheeler 2024) to
find the dominant tandem repeat, handling substitutions and small indels
that a naive k-mer scan would miss. It emits one row per intron that
contains at least one repeat passing the configured cutoffs:

  dominant_consensus       ULTRA's estimated repeat consensus (may be a rotation)
  dominant_canonical       lex-min rotation of consensus and its RC — collapses
                           all rotations + both strands to one identifier
  period                   repeat unit length (bp)
  score                    ULTRA HMM score
  copies                   estimated repeat unit count
  substitutions            mismatches vs the perfect array (degeneracy signal)
  insertions               small insertions
  deletions                small deletions
  covered_bp               End - Start of the dominant repeat
  cover_frac               covered_bp / intron_len
  telomere_match           True if canonical matches a config telomere motif
  telomere_match_name      the matching config motif (e.g. TTAGGG), else ""

Intended for MANUAL CURATION — no filtering by "known telomere" is done here.

Wildcarded per-genome mode: --single-genome <gid> processes just that
genome, so the Snakefile can fan out one sbatch job per genome.
"""
from __future__ import annotations
import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import rc, rotations, read_fasta, find_genome_files, open_maybe_gz

# scan_telotrons contains the gt-driven intron extractor; reuse it verbatim
# so both scans agree on which introns exist for a given GFF (silent
# gt-truncation trap already handled there).
from scan_telotrons import run_gt_introns, load_introns


def canonical_unit(unit: str) -> str:
    """Lex-min rotation of unit and of its reverse complement.

    Collapses every rotation and both strands of the same tandem unit to a
    single identifier — TTAGGG, TAGGGT, GGGTTA, CCCTAA all become AACCCT.
    """
    u = unit.upper()
    if not u:
        return ""
    both = rotations(u) + rotations(rc(u))
    return min(both)


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


INTRON_HEADER_SEP = "|"


def _encode_intron_id(idx: int, r) -> str:
    # Encode the intron row index inside the FASTA header so we can join
    # ULTRA's per-record output back to the intron table without carrying a
    # separate side-file. Row index is enough — everything else lives in
    # `introns`.
    return f"i{idx}"


def _run_ultra(fa_path: str, out_tsv: str, opts: dict, tmp_dir: str,
               ultra_bin: str, threads: int) -> None:
    """Run ULTRA and write its TSV output to out_tsv. Raises on non-zero exit."""
    # ULTRA streams to stdout in --tsv mode when -o is unset.
    cmd = [
        ultra_bin, "--tsv", "-c",                # counts columns: subs/ins/del
        "--min_unit", str(opts["min_units"]),
        "--min_length", str(opts["min_length"]),
        "-p", str(opts["max_period"]),
        "-s", str(opts["min_score"]),
        "-t", str(threads),
        fa_path,
    ]
    with open(out_tsv, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or "")[-500:].replace("\n", " | ")
        raise RuntimeError(f"ULTRA exit {proc.returncode}: {tail}")


def _parse_ultra_tsv(tsv_path: str) -> dict[str, list[dict]]:
    """Group ULTRA rows by SeqID; keep only the numeric columns we consume.

    ULTRA's TSV header (v1.2.1): SeqID, Start, End, Period, Score, Consensus,
    #copies, #substitutions, #insertions, #deletions, #Subrepeats,
    SubrepeatStarts, SubrepeatConsensi, [Sequence].
    """
    by_seq: dict[str, list[dict]] = {}
    with open(tsv_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("=") or line.startswith("*"):
                continue
            if line.startswith("SeqID"):
                header = line.split("\t")
                continue
            if header is None:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            row = dict(zip(header, parts))
            try:
                rec = {
                    "start": int(row["Start"]),
                    "end": int(row["End"]),
                    "period": int(row["Period"]),
                    "score": float(row["Score"]),
                    "consensus": row.get("Consensus", "").upper(),
                    "copies": float(row.get("#copies", "0") or 0),
                    "subs": int(float(row.get("#substitutions", "0") or 0)),
                    "ins": int(float(row.get("#insertions", "0") or 0)),
                    "dels": int(float(row.get("#deletions", "0") or 0)),
                }
            except (KeyError, ValueError):
                continue
            by_seq.setdefault(row["SeqID"], []).append(rec)
    return by_seq


COLUMNS = [
    "genome_id", "organism", "group", "source",
    "seqid", "start", "end", "strand", "tx_id", "gene_id", "intron_index", "intron_len",
    "dominant_consensus", "dominant_canonical", "period",
    "score", "copies", "substitutions", "insertions", "deletions",
    "covered_bp", "cover_frac",
    "telomere_match", "telomere_match_name",
]


def scan_one(row, opts, ultra_bin: str, threads: int):
    gid, org, group, source = row
    try:
        fa, gff = find_genome_files(gid, opts["refseq_dir"], opts["tara_dir"])
    except FileNotFoundError:
        return {"error": f"missing fasta or gff for {gid}", "rows": []}
    seqs = read_fasta(fa)

    with tempfile.TemporaryDirectory(prefix=f"repeat_{gid}_") as tmp:
        # 1. Build intron table via the same gt pipeline scan_telotrons uses.
        gff_with_introns = os.path.join(tmp, "with_introns.gff3")
        try:
            run_gt_introns(gff, gff_with_introns)
        except Exception as e:
            return {"error": f"gt failed for {gid}: {e}", "rows": []}
        introns = load_introns(gff_with_introns)
        if introns.empty:
            return {"error": "", "rows": []}
        introns = introns.reset_index(drop=True)

        # 2. Emit an intron FASTA — one record per intron, strand-corrected so
        #    ULTRA's Consensus is reported on the transcript strand.
        intron_fa = os.path.join(tmp, "introns.fa")
        keep: dict[str, tuple] = {}
        with open(intron_fa, "w") as fh:
            for idx, r in enumerate(introns.itertuples(index=False)):
                contig = seqs.get(r.seqid)
                if contig is None:
                    continue
                s0 = int(r.start) - 1
                e0 = int(r.end)
                if e0 > len(contig) or e0 - s0 < opts["min_intron_len"]:
                    continue
                raw = contig[s0:e0]
                seq = revcomp_if_neg(raw, r.strand)
                rid = _encode_intron_id(idx, r)
                keep[rid] = (idx, r, len(seq))
                fh.write(f">{rid}\n{seq}\n")
        if not keep:
            return {"error": "", "rows": []}

        # 3. Run ULTRA on the intron FASTA.
        ultra_tsv = os.path.join(tmp, "ultra.tsv")
        try:
            _run_ultra(intron_fa, ultra_tsv, opts, tmp, ultra_bin, threads)
        except RuntimeError as e:
            return {"error": f"{gid}: {e}", "rows": []}

        # 4. For each intron, pick the ULTRA record with the largest covered_bp
        #    (ties broken by score). ULTRA can report multiple repeats per
        #    sequence (subrepeat structure); the largest wins.
        by_seq = _parse_ultra_tsv(ultra_tsv)

    telomere_lookup = opts["telomere_lookup"]
    min_frac = opts["min_frac"]

    rows = []
    for rid, (idx, r, ilen) in keep.items():
        records = by_seq.get(rid) or []
        if not records:
            continue
        rec = max(records, key=lambda x: (x["end"] - x["start"], x["score"]))
        covered = rec["end"] - rec["start"]
        frac = covered / ilen if ilen else 0.0
        if frac < min_frac:
            continue
        canon = canonical_unit(rec["consensus"])
        match_name = telomere_lookup.get(canon, "")
        rows.append([
            gid, org, group, source,
            r.seqid, int(r.start), int(r.end), r.strand,
            r.tx_id, r.gene_id, r.intron_index, ilen,
            rec["consensus"], canon, rec["period"],
            round(rec["score"], 3), round(rec["copies"], 3),
            rec["subs"], rec["ins"], rec["dels"],
            covered, round(frac, 4),
            bool(match_name), match_name,
        ])
    return {"error": "", "rows": rows}


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
    ap.add_argument("--ultra-bin", default="ultra",
                    help="ULTRA binary path/name (default: on PATH)")
    ap.add_argument("--min-length", type=int, default=20,
                    help="ULTRA --min_length: minimum reportable repeat length (bp)")
    ap.add_argument("--min-units", type=int, default=3,
                    help="ULTRA --min_unit: minimum reportable repeat unit count")
    ap.add_argument("--max-period", type=int, default=20,
                    help="ULTRA -p: maximum detectable repeat period (bp). "
                         "20 comfortably covers telomere-scale repeats (4-8 bp) "
                         "and their small tandem oligomer variants.")
    ap.add_argument("--min-score", type=float, default=0.0,
                    help="ULTRA -s: minimum HMM score (default 0 = permissive; "
                         "curate downstream)")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="post-ULTRA filter: only emit introns whose dominant "
                         "repeat covers this fraction of the intron")
    ap.add_argument("--min-intron-len", type=int, default=30)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    if not shutil.which(args.ultra_bin):
        sys.exit(f"ULTRA binary not on PATH: {args.ultra_bin!r} — "
                 f"conda install -c bioconda ultra")

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
        "min_length": args.min_length,
        "min_units": args.min_units,
        "max_period": args.max_period,
        "min_score": args.min_score,
        "min_frac": args.min_frac,
        "min_intron_len": args.min_intron_len,
    }

    # Per-genome: single ULTRA process, ULTRA itself parallelizes across cores.
    n_err = 0
    n_rows = 0
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(COLUMNS)
        for row in rows:
            res = scan_one(row, opts, args.ultra_bin, args.threads)
            if res["error"]:
                n_err += 1
                print(f"[repeat_scan] {res['error']}", file=sys.stderr)
            for r in res["rows"]:
                w.writerow(r)
                n_rows += 1
    print(f"[repeat_scan] wrote {n_rows} rows across {len(rows)} genomes "
          f"({n_err} errors) → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
