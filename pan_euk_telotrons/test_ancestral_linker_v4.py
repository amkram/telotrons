#!/usr/bin/env python3
"""
v4: extend the ancestral-linker test to ALL fill_in loci with >=3 cross-species
orthologs. Reads sequences directly from
linker_analysis/orthology_v2/results/ortholog_results.json (no need to rerun
the alignment pipeline).

For each host telotron locus:
  - Strip telomeric kmers from host intron -> linker
  - For each cross-species target, BLAST host_linker vs target.intron_dna
  - Tally hits, report aggregate per-locus and global

Lineage-aware kmer set (TTAGGG vertebrate vs TTTAGGG plant/protozoan).
"""

import json
import subprocess
import tempfile
import csv
from pathlib import Path
from collections import Counter, defaultdict

ORTHO_JSON = Path("/scratch1/alex/telotrons/pan_euk_telotrons/downstream/linker_analysis/orthology_v2/results/ortholog_results.json")
VERDICT_TSV = Path("/scratch1/alex/telotrons/pan_euk_telotrons/downstream/linker_analysis/orthology_v2/results/locus_verdict_summary.tsv")
OUT_DIR = Path("/scratch1/alex/telotrons/pan_euk_telotrons/downstream/linker_analysis/orthology_v2/results")

# Both 6mer (vertebrate/fungal) and 7mer (plant/protozoan) telomeric kmer sets.
# We'll mask anything matching either family on either strand.
TELO = set()
for unit in ["TTAGGG", "TTTAGGG"]:
    rotations = [unit[i:] + unit[:i] for i in range(len(unit))]
    TELO.update(rotations)
    rev = unit[::-1].translate(str.maketrans("ACGT", "TGCA"))
    rotations = [rev[i:] + rev[:i] for i in range(len(rev))]
    TELO.update(rotations)


def telo_mask(s):
    s = s.upper()
    n = len(s)
    m = [0] * n
    i = 0
    while i < n:
        hit = False
        for k in TELO:
            if s[i:i+len(k)] == k:
                for j in range(i, min(i+len(k), n)):
                    m[j] = 1
                i += len(k); hit = True
                break
        if not hit:
            i += 1
    return m


def low_complexity_frac(seq, k=4):
    """Fraction of bases covered by the most common k-mer + its 1mm neighbors.
    High value = microsatellite/repetitive. Returns 0.0 for normal sequence."""
    if len(seq) < 2 * k:
        return 0.0
    from collections import Counter
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    if not kmers:
        return 0.0
    c = Counter(kmers)
    top = c.most_common(1)[0][0]
    # Build 1-mismatch neighborhood of top
    nbr = {top}
    for i in range(k):
        for b in "ACGT":
            if b != top[i]:
                nbr.add(top[:i] + b + top[i+1:])
    n_in_nbr = sum(c[k] for k in nbr if k in c)
    return n_in_nbr / max(1, len(kmers))


def is_telomicrosatellite(seq):
    """Return True if seq is dominated by telo-family repeats with 1mm wobble
    (i.e. the kind of degenerate-telomeric content that masquerades as 'linker'
    after exact-kmer stripping)."""
    if len(seq) < 30:
        return False
    return low_complexity_frac(seq, k=4) >= 0.50


def blast(q, s):
    if len(q) < 15 or len(s) < 15:
        return []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td/"q.fa").write_text(f">q\n{q}\n")
        (td/"s.fa").write_text(f">s\n{s}\n")
        try:
            r = subprocess.run(
                ["blastn", "-task", "blastn", "-query", str(td/"q.fa"),
                 "-subject", str(td/"s.fa"),
                 "-evalue", "1", "-word_size", "7", "-dust", "no",
                 "-outfmt", "6 pident length evalue qstart qend sstart send"],
                capture_output=True, text=True, timeout=20,
            )
        except Exception:
            return []
    hits = []
    for line in r.stdout.strip().split("\n"):
        if not line:
            continue
        p = line.split("\t")
        hits.append({"pident": float(p[0]), "length": int(p[1]),
                     "evalue": float(p[2])})
    return hits


def best_hit(hits, min_pid=70, min_len=20):
    qual = [h for h in hits if h["pident"] >= min_pid and h["length"] >= min_len]
    if not qual:
        return None
    return max(qual, key=lambda x: x["pident"] * x["length"])


def main():
    # Load fill_in locus IDs
    fillin = set()
    with open(VERDICT_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["bin"] == "fill_in":
                fillin.add(r["locus"])

    # Load ortholog data
    data = json.load(open(ORTHO_JSON))

    # Match locus_id -> entry
    entries = []
    n_excluded_1mm = 0
    for r in data:
        locus_id = f"{r['host_sp'].split('__')[0]}__{r['host_mrna']}__{r['host_start']}_{r['host_end']}"
        if locus_id not in fillin:
            continue
        targets = r.get("targets", {})
        if len(targets) < 3:
            continue
        # Exclude 1-mismatch ULTRA hits — CLAUDE.md flags these as ~10:1 noise
        telo_type = r.get("telo_type", "")
        if "~1mm" in telo_type or "1mm" in telo_type:
            n_excluded_1mm += 1
            continue
        entries.append((locus_id, r))

    print(f"fill_in loci with >=3 targets (exact telo only): {len(entries)}")
    print(f"  Excluded {n_excluded_1mm} loci with 1-mismatch telo_type (likely microsatellites, not telotrons)")
    # n_excluded_lowcomplex reported below after main loop

    rows = []
    locus_summary = []
    n_excluded_lowcomplex = 0

    for locus_id, r in entries:
        host_intron = r["host_intron_seq"].upper()
        m = telo_mask(host_intron)
        host_linker = "".join(c for c, mb in zip(host_intron, m) if mb == 0)
        host_telo_frac = sum(m) / len(m) if m else 0
        # Use k=4 AND k=5 dominance — k=5 catches CTAAC/CCTAA degenerate-telo
        # microsatellites that survive the exact-kmer strip; k=4 catches GGGAT
        # microsatellites that survive 1mm-telo filtering.
        host_lc4 = low_complexity_frac(host_linker, k=4) if len(host_linker) >= 30 else 0.0
        host_lc5 = low_complexity_frac(host_linker, k=5) if len(host_linker) >= 30 else 0.0
        if host_lc4 >= 0.40 or host_lc5 >= 0.40:
            n_excluded_lowcomplex += 1
            continue

        n_with_strong_blast = 0
        n_plain_targets = 0  # targets with low telo content (clean test)
        per_target = []

        for target_id, tinfo in r["targets"].items():
            t_intron = (tinfo.get("intron_dna") or "").upper()
            if not t_intron:
                continue
            tm = telo_mask(t_intron)
            t_telo_frac = sum(tm) / len(tm) if tm else 0
            t_linker = "".join(c for c, mb in zip(t_intron, tm) if mb == 0)

            hits_full = blast(host_linker, t_intron)
            best_full = best_hit(hits_full, min_pid=70, min_len=20)

            row = {
                "locus": locus_id,
                "host": r["host_sp"],
                "host_intron_len": len(host_intron),
                "host_telo_frac": round(host_telo_frac, 3),
                "host_linker_len": len(host_linker),
                "target": target_id,
                "target_intron_len": len(t_intron),
                "target_telo_frac": round(t_telo_frac, 3),
                "target_linker_len": len(t_linker),
                "blast_best_pid": round(best_full["pident"], 1) if best_full else None,
                "blast_best_len": best_full["length"] if best_full else None,
                "blast_n_hits": len(hits_full),
            }
            rows.append(row)
            per_target.append(row)
            if t_telo_frac < 0.10 and len(host_linker) >= 25:
                n_plain_targets += 1
                if best_full:
                    n_with_strong_blast += 1

        locus_summary.append({
            "locus": locus_id,
            "host": r["host_sp"],
            "n_targets": len(r["targets"]),
            "host_intron_len": len(host_intron),
            "host_linker_len": len(host_linker),
            "host_telo_frac": round(host_telo_frac, 3),
            "n_plain_targets": n_plain_targets,
            "n_with_strong_blast": n_with_strong_blast,
            "verdict": "ANCESTRAL_HOMOLOGY" if n_with_strong_blast >= 2
                      else ("WEAK" if n_with_strong_blast == 1 else "NO_HOMOLOGY"),
        })

    print(f"  Excluded {n_excluded_lowcomplex} loci where host linker is microsatellite (lowcomp>=0.50)")
    print(f"  Net loci tested: {len(locus_summary)}")

    # Sort summary: ANCESTRAL first
    rank = {"ANCESTRAL_HOMOLOGY": 0, "WEAK": 1, "NO_HOMOLOGY": 2}
    locus_summary.sort(key=lambda x: (rank[x["verdict"]], -x["n_with_strong_blast"]))

    # Write outputs
    cols = list(rows[0].keys())
    with open(OUT_DIR / "ancestral_linker_v4.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) if r[c] is not None else "" for c in cols) + "\n")

    cols2 = list(locus_summary[0].keys())
    with open(OUT_DIR / "ancestral_linker_v4_per_locus.tsv", "w") as f:
        f.write("\t".join(cols2) + "\n")
        for r in locus_summary:
            f.write("\t".join(str(r[c]) for c in cols2) + "\n")

    # Aggregate
    n_total = len(locus_summary)
    n_anc = sum(1 for r in locus_summary if r["verdict"] == "ANCESTRAL_HOMOLOGY")
    n_weak = sum(1 for r in locus_summary if r["verdict"] == "WEAK")
    n_none = sum(1 for r in locus_summary if r["verdict"] == "NO_HOMOLOGY")

    print(f"\n=== Per-locus verdicts ===")
    print(f"  ANCESTRAL_HOMOLOGY (≥2 plain targets with strong BLAST): {n_anc}/{n_total}")
    print(f"  WEAK (exactly 1 plain target with strong BLAST):         {n_weak}/{n_total}")
    print(f"  NO_HOMOLOGY (0 plain targets with strong BLAST):         {n_none}/{n_total}")

    print(f"\n--- ANCESTRAL_HOMOLOGY loci ---")
    for r in locus_summary:
        if r["verdict"] != "ANCESTRAL_HOMOLOGY":
            continue
        print(f"  {r['host'][:35]:>35s}  intron={r['host_intron_len']:>5d}bp  "
              f"linker={r['host_linker_len']:>4d}bp  "
              f"plain_n={r['n_plain_targets']}  strong_n={r['n_with_strong_blast']}  "
              f"{r['locus']}")

    print(f"\n--- WEAK loci ---")
    for r in locus_summary:
        if r["verdict"] != "WEAK":
            continue
        print(f"  {r['host'][:35]:>35s}  intron={r['host_intron_len']:>5d}bp  "
              f"linker={r['host_linker_len']:>4d}bp  "
              f"plain_n={r['n_plain_targets']}  strong_n={r['n_with_strong_blast']}")

    print(f"\n--- NO_HOMOLOGY loci (showing first 15) ---")
    for r in locus_summary:
        if r["verdict"] != "NO_HOMOLOGY":
            continue
        print(f"  {r['host'][:35]:>35s}  intron={r['host_intron_len']:>5d}bp  "
              f"linker={r['host_linker_len']:>4d}bp  "
              f"plain_n={r['n_plain_targets']}")

    # Save summary JSON
    summary = {
        "n_loci": n_total,
        "n_ancestral_homology": n_anc,
        "n_weak": n_weak,
        "n_no_homology": n_none,
        "max_blast_length_observed": max((r["blast_best_len"] for r in rows if r["blast_best_len"]), default=0),
        "max_blast_pident_observed": max((r["blast_best_pid"] for r in rows if r["blast_best_pid"]), default=0),
        "ancestral_loci": [r["locus"] for r in locus_summary if r["verdict"] == "ANCESTRAL_HOMOLOGY"],
        "weak_loci": [r["locus"] for r in locus_summary if r["verdict"] == "WEAK"],
    }
    json.dump(summary, open(OUT_DIR / "ancestral_linker_v4_summary.json", "w"), indent=2)
    print(f"\n=== Summary JSON ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
