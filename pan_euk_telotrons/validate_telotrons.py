#!/usr/bin/env python3
"""
Validate and filter candidate telotrons from pan-eukaryotic survey.

Checks applied:
  1. Purity ≥85% telomeric hexamer coverage
  2. GT-AG (or CT-AC reverse complement) canonical splice sites
  3. Minimum length ≥30bp
  4. Strand consistency (template vs coding vs converging)
  5. Boundary motif analysis (TERC fingerprint: CTAA/AACC enrichment)
  6. Repeat unit identification (TTAGGG vs CCCTAA vs variants)
  7. Clustering — do telotrons cluster on specific contigs?
  8. Array purity profile — truly tandem vs scattered hexamers
  9. Splice signal intrinsic to repeat — can the repeat provide GT/AG?
 10. Flag potential assembly artifacts (identical sequences, tiny contigs)

Usage:
  python3 validate_telotrons.py                  # process all TSV files in telotrons/
  python3 validate_telotrons.py telotrons/GCF_001643675.1_*.tsv  # specific file(s)
"""

import csv, json, os, re, sys
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Telomeric hexamer sets
# ---------------------------------------------------------------------------
def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

TELO_HEX_6 = _rotations("TTAGGG") | _rotations("CCCTAA")

FWD_PATS = sorted(_rotations("TTAGGG"))  # TTAGGG rotations
REV_PATS = sorted(_rotations("CCCTAA"))  # CCCTAA rotations


def rc(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


# ---------------------------------------------------------------------------
# Validation functions
# ---------------------------------------------------------------------------

def check_tandem_structure(seq, hexset=TELO_HEX_6):
    """Check if the telomeric content forms a contiguous tandem array
    vs scattered hexamers. Returns (longest_run_bp, longest_run_frac, n_gaps)."""
    n = len(seq)
    if n < 6:
        return 0, 0.0, 0

    covered = bytearray(n)
    for h in hexset:
        start = 0
        while True:
            idx = seq.find(h, start)
            if idx == -1:
                break
            for j in range(idx, min(idx + len(h), n)):
                covered[j] = 1
            start = idx + 1

    # Find longest contiguous run
    runs = []
    run_start = None
    for i in range(n):
        if covered[i]:
            if run_start is None:
                run_start = i
        else:
            if run_start is not None:
                runs.append((run_start, i))
                run_start = None
    if run_start is not None:
        runs.append((run_start, n))

    if not runs:
        return 0, 0.0, 0

    longest = max(runs, key=lambda x: x[1] - x[0])
    longest_len = longest[1] - longest[0]
    n_gaps = len(runs) - 1  # gaps between runs

    return longest_len, longest_len / n, n_gaps


def check_splice_from_repeat(seq, strand_orient):
    """Check if GT donor and AG acceptor can be derived from the telomeric
    repeat itself (the key telotron property).

    TTAGGG contains AG at positions 2-3 (acceptor) but no GT.
    GT can come from cross-boundary: ...GGG|TTA... → GT at junction not standard.

    For the insertion to function as an intron:
    - Need GT at 5' end and AG at 3' end
    - Or CT at 5' and AC at 3' (reverse complement on - strand)
    """
    d2 = seq[:2].upper()
    a2 = seq[-2:].upper()

    # Normalize to sense strand
    if d2 == "CT" and a2 == "AC":
        norm_seq = rc(seq)
        d2 = norm_seq[:2]
        a2 = norm_seq[-2:]
    else:
        norm_seq = seq

    result = {
        "donor": d2,
        "acceptor": a2,
        "canonical": d2 == "GT" and a2 == "AG",
        "gt_from_repeat": False,
        "ag_from_repeat": False,
    }

    # Check if GT comes from TTAGGG boundary: ...TTAGG|GTTAG...
    # The canonical TTAGGG repeat has GT at rotation GTTAGG
    if d2 == "GT":
        # GT at position 0 = rotation GTTAGG = within repeat
        result["gt_from_repeat"] = norm_seq[:6].upper() in TELO_HEX_6

    # AG at end: TTAGGG has AG at positions 2-3
    if a2 == "AG":
        result["ag_from_repeat"] = norm_seq[-6:].upper() in TELO_HEX_6

    return result


def check_boundary_motifs(seq, orient):
    """Extract boundary tetranucleotides for TERC fingerprint analysis."""
    norm = seq.upper()
    d2 = norm[:2]
    a2 = norm[-2:]

    # Normalize to GT-AG
    if d2 == "CT" and a2 == "AC":
        norm = rc(norm)
        d2 = norm[:2]
        a2 = norm[-2:]

    result = {"donor_4": None, "acceptor_4": None}
    if len(norm) >= 8 and d2 == "GT" and a2 == "AG":
        result["donor_4"] = norm[2:6]
        result["acceptor_4"] = norm[-6:-2]
    return result


def check_artifact_signals(seq, contig, start, end, all_seqs):
    """Flag potential assembly artifacts."""
    flags = []

    # Exact duplicate sequences (different loci, same sequence)
    seq_upper = seq.upper()
    count = all_seqs.get(seq_upper, 0)
    if count > 1:
        flags.append(f"DUPLICATE_SEQ({count}x)")

    # Very short contig (could be unplaced fragment)
    # We don't have contig lengths here, but flag if coords suggest tiny contig
    if end < 5000 and start < 1000:
        flags.append("NEAR_CONTIG_START")

    # Suspiciously regular — perfectly repeating with no variation
    if len(seq) >= 12:
        period6 = seq_upper[:6]
        perfect_repeats = sum(1 for i in range(0, len(seq_upper) - 5, 6)
                             if seq_upper[i:i+6] == period6)
        total_possible = len(seq_upper) // 6
        if total_possible > 0 and perfect_repeats / total_possible > 0.95:
            flags.append("PERFECTLY_PERIODIC")

    return flags


def classify_repeat_unit(seq):
    """Determine which telomeric repeat unit dominates."""
    seq = seq.upper()
    counts = {}
    for pat in FWD_PATS:
        counts[pat] = len(re.findall(pat, seq))
    for pat in REV_PATS:
        counts[pat] = len(re.findall(pat, seq))

    fwd_total = sum(counts.get(p, 0) for p in FWD_PATS)
    rev_total = sum(counts.get(p, 0) for p in REV_PATS)

    if fwd_total + rev_total == 0:
        return "unknown", 0, 0

    if fwd_total > rev_total:
        dominant = "TTAGGG"
        top_rot = max(FWD_PATS, key=lambda p: counts[p])
    else:
        dominant = "CCCTAA"
        top_rot = max(REV_PATS, key=lambda p: counts[p])

    return dominant, fwd_total, rev_total


# ---------------------------------------------------------------------------
# Main validation pipeline
# ---------------------------------------------------------------------------

def validate_file(tsv_path):
    """Validate all candidate telotrons in a TSV file."""
    species = Path(tsv_path).stem.replace("_telotrons", "")

    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    if not rows:
        return None

    # Count duplicate sequences for artifact detection
    all_seqs = Counter()
    for r in rows:
        seq = r.get("intron_seq", "").upper()
        if seq:
            all_seqs[seq] += 1

    validated = []
    rejected = []

    for r in rows:
        seq = r.get("intron_seq", "")
        if not seq:
            continue

        cov = float(r.get("telo_coverage", 0))
        length = int(r.get("length", 0))
        donor = r.get("donor", "")
        acceptor = r.get("acceptor", "")
        orient = r.get("orientation", "")
        contig = r.get("contig", "")
        start = int(r.get("start", 0))
        end = int(r.get("end", 0))

        record = {
            "species": species,
            "contig": contig,
            "start": start,
            "end": end,
            "length": length,
            "strand": r.get("strand", ""),
            "gene": r.get("mrna", ""),
            "donor": donor,
            "acceptor": acceptor,
            "telo_coverage": round(cov, 4),
            "orientation": orient,
            "fwd_frac": float(r.get("fwd_frac", 0.5)),
        }

        # --- Apply filters ---
        reject_reasons = []

        # 1. Purity
        if cov < 0.85:
            reject_reasons.append(f"LOW_PURITY({cov:.1%})")

        # 2. Splice sites
        is_canonical = (donor == "GT" and acceptor == "AG") or \
                       (donor == "CT" and acceptor == "AC")
        if not is_canonical:
            reject_reasons.append(f"NON_CANONICAL_SPLICE({donor}-{acceptor})")

        # 3. Minimum length
        if length < 30:
            reject_reasons.append(f"TOO_SHORT({length}bp)")

        # 4. Tandem structure
        longest_run, run_frac, n_gaps = check_tandem_structure(seq)
        record["longest_telo_run"] = longest_run
        record["telo_run_frac"] = round(run_frac, 3)
        record["n_gaps"] = n_gaps
        if run_frac < 0.50:
            reject_reasons.append(f"SCATTERED_HEXAMERS(run={run_frac:.0%})")

        # 5. Splice signals from repeat
        splice_info = check_splice_from_repeat(seq, orient)
        record["splice_canonical"] = splice_info["canonical"]
        record["gt_from_repeat"] = splice_info["gt_from_repeat"]
        record["ag_from_repeat"] = splice_info["ag_from_repeat"]

        # 6. Boundary motifs
        bnd = check_boundary_motifs(seq, orient)
        record["donor_4"] = bnd["donor_4"]
        record["acceptor_4"] = bnd["acceptor_4"]

        # 7. Repeat unit
        dominant, fwd_n, rev_n = classify_repeat_unit(seq)
        record["repeat_unit"] = dominant
        record["fwd_hexamers"] = fwd_n
        record["rev_hexamers"] = rev_n

        # 8. Artifact checks
        flags = check_artifact_signals(seq, contig, start, end, all_seqs)
        record["artifact_flags"] = flags

        # 9. Sequence (truncated for output)
        record["seq_preview"] = seq[:80]

        if reject_reasons:
            record["reject_reasons"] = reject_reasons
            rejected.append(record)
        else:
            record["reject_reasons"] = []
            validated.append(record)

    return {
        "species": species,
        "file": str(tsv_path),
        "total_candidates": len(rows),
        "validated": len(validated),
        "rejected": len(rejected),
        "validated_records": validated,
        "rejected_records": rejected,
    }


# ---------------------------------------------------------------------------
# Cross-species summary
# ---------------------------------------------------------------------------

def write_summary(all_results, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect all validated telotrons
    all_validated = []
    species_summary = []

    for result in all_results:
        if result is None:
            continue
        v = result["validated_records"]
        all_validated.extend(v)

        if result["validated"] > 0:
            # Per-species stats
            orientations = Counter(r["orientation"] for r in v)
            splice_from_repeat = sum(1 for r in v if r["gt_from_repeat"] and r["ag_from_repeat"])
            donor_4s = Counter(r["donor_4"] for r in v if r["donor_4"])
            acceptor_4s = Counter(r["acceptor_4"] for r in v if r["acceptor_4"])
            contigs = Counter(r["contig"] for r in v)

            species_summary.append({
                "species": result["species"],
                "n_validated": result["validated"],
                "n_rejected": result["rejected"],
                "n_total": result["total_candidates"],
                "orientations": dict(orientations),
                "splice_from_repeat": splice_from_repeat,
                "top_donor_4": donor_4s.most_common(3),
                "top_acceptor_4": acceptor_4s.most_common(3),
                "n_contigs": len(contigs),
                "top_contigs": contigs.most_common(5),
                "mean_length": round(sum(r["length"] for r in v) / len(v), 1) if v else 0,
                "mean_purity": round(sum(r["telo_coverage"] for r in v) / len(v), 4) if v else 0,
                "artifact_flags": Counter(f for r in v for f in r["artifact_flags"]),
            })

    species_summary.sort(key=lambda x: -x["n_validated"])

    # --- Markdown report ---
    md_path = output_dir / "validated_telotrons_report.md"
    with open(md_path, 'w') as f:
        f.write("# Validated Telotrons — Pan-Eukaryotic Survey\n\n")
        f.write(f"**Date:** {datetime.now().isoformat()[:10]}\n")
        f.write(f"**Total validated telotrons:** {len(all_validated)}\n")
        f.write(f"**Species with validated telotrons:** {len(species_summary)}\n\n")

        f.write("## Validation Criteria\n\n")
        f.write("1. Telomeric coverage ≥85%\n")
        f.write("2. Canonical splice sites (GT-AG or CT-AC)\n")
        f.write("3. Length ≥30bp\n")
        f.write("4. Contiguous tandem array (longest run ≥50% of intron)\n\n")

        f.write("## Species Summary\n\n")
        f.write("| Species | Validated | Rejected | Mean length | Mean purity | "
                "Template | Coding | Converging | Top 3' motif | Clustered? |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---|---|\n")

        for s in species_summary:
            o = s["orientations"]
            top_a = s["top_acceptor_4"][0] if s["top_acceptor_4"] else ("—", 0)
            top_a_str = f"{top_a[0]}({top_a[1]})" if isinstance(top_a, tuple) else str(top_a)
            clustered = "YES" if s["n_contigs"] < s["n_validated"] / 3 else "no"
            f.write(
                f"| {s['species'][:40]} "
                f"| {s['n_validated']} "
                f"| {s['n_rejected']} "
                f"| {s['mean_length']:.0f}bp "
                f"| {s['mean_purity']:.1%} "
                f"| {o.get('template-strand', 0)} "
                f"| {o.get('coding-strand', 0)} "
                f"| {o.get('converging', 0)} "
                f"| {top_a_str} "
                f"| {clustered} |\n"
            )

        # Detailed per-species sections
        for s in species_summary:
            f.write(f"\n### {s['species']}\n\n")
            f.write(f"- **Validated:** {s['n_validated']} / {s['n_total']}\n")
            f.write(f"- **Mean length:** {s['mean_length']:.0f} bp\n")
            f.write(f"- **Mean purity:** {s['mean_purity']:.1%}\n")
            f.write(f"- **Orientations:** {dict(s['orientations'])}\n")
            f.write(f"- **Splice signals from repeat:** {s['splice_from_repeat']}/{s['n_validated']}\n")
            f.write(f"- **Top 3' boundary motifs:** {s['top_acceptor_4']}\n")
            f.write(f"- **Top 5' boundary motifs:** {s['top_donor_4']}\n")
            f.write(f"- **Contigs with telotrons:** {s['n_contigs']} "
                    f"(top: {s['top_contigs'][:3]})\n")
            if s["artifact_flags"]:
                f.write(f"- **Artifact flags:** {dict(s['artifact_flags'])}\n")
            f.write("\n")

    # --- JSON with all validated records ---
    json_path = output_dir / "validated_telotrons.json"
    with open(json_path, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "total_validated": len(all_validated),
            "species_summary": species_summary,
            "validated_telotrons": all_validated[:2000],  # cap for file size
        }, f, indent=2, default=str)

    # --- TSV of validated telotrons ---
    tsv_path = output_dir / "validated_telotrons.tsv"
    cols = ["species", "contig", "start", "end", "length", "strand", "gene",
            "donor", "acceptor", "telo_coverage", "orientation", "repeat_unit",
            "longest_telo_run", "telo_run_frac", "n_gaps",
            "splice_canonical", "gt_from_repeat", "ag_from_repeat",
            "donor_4", "acceptor_4", "artifact_flags", "seq_preview"]
    with open(tsv_path, 'w') as f:
        f.write("\t".join(cols) + "\n")
        for r in all_validated:
            vals = [str(r.get(c, "")) for c in cols]
            f.write("\t".join(vals) + "\n")

    print(f"\nOutputs:")
    print(f"  {md_path}")
    print(f"  {json_path}")
    print(f"  {tsv_path}")

    # Print summary to stdout
    print(f"\n{'='*70}")
    print(f"VALIDATED TELOTRONS SUMMARY")
    print(f"{'='*70}")
    print(f"Total validated: {len(all_validated)}")
    print(f"Species with validated telotrons: {len(species_summary)}")
    print()
    for s in species_summary:
        o = s["orientations"]
        print(f"  {s['species']:<40s} {s['n_validated']:>4d} telotrons  "
              f"(T:{o.get('template-strand',0)} C:{o.get('coding-strand',0)} "
              f"X:{o.get('converging',0)})  "
              f"mean={s['mean_length']:.0f}bp  "
              f"purity={s['mean_purity']:.1%}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    telotron_dir = Path(__file__).parent / "telotrons"
    results_dir = Path(__file__).parent / "validated_results"

    if len(sys.argv) > 1:
        tsv_files = [Path(f) for f in sys.argv[1:]]
    else:
        tsv_files = sorted(telotron_dir.glob("*_telotrons.tsv"))

    if not tsv_files:
        print("No telotron TSV files found.")
        sys.exit(1)

    print(f"Validating {len(tsv_files)} telotron files...\n")

    all_results = []
    for tsv in tsv_files:
        result = validate_file(tsv)
        if result and result["total_candidates"] > 0:
            v = result["validated"]
            r = result["rejected"]
            t = result["total_candidates"]
            species = result["species"][:45]
            if v > 0:
                print(f"  ★ {species:<45s} {v:>4d}/{t} validated")
            else:
                print(f"    {species:<45s} {r:>4d}/{t} rejected")
            all_results.append(result)

    write_summary(all_results, results_dir)


if __name__ == "__main__":
    main()
