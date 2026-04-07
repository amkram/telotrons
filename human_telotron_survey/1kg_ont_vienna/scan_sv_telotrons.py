#!/usr/bin/env python3
"""
Scan 1KG ONT Vienna SV data for telomeric-repeat insertions.

Two modes:
  1. VCF mode: scan inserted ALT sequences in the merged VCF for TTAGGG arrays
  2. SVTig mode: scan SVTig FASTA contigs for TTAGGG arrays

Reports any SV insertion with ≥50bp of ≥85% telomeric-repeat content,
plus a broader scan at lower thresholds for context.
"""

import gzip, re, sys, json, os
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# Telomeric hexamer sets
# ---------------------------------------------------------------------------
def _rotations(base):
    return {base[i:] + base[:i] for i in range(len(base))}

TELO_HEX = _rotations("TTAGGG") | _rotations("CCCTAA")

FWD_PATS = ["TTAGGG", "TAGGGT", "AGGGTT", "GGGTTA", "GGTTAG", "GTTAGG"]
REV_PATS = ["CCCTAA", "CCTAAC", "CTAACC", "TAACCC", "AACCCT", "ACCCTA"]


def telo_coverage(seq, hexset=TELO_HEX):
    """Fraction of bases covered by any hexamer in set."""
    n = len(seq)
    if n < 6:
        return 0.0
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
    return sum(covered) / n


def strand_ratio(seq):
    fwd = sum(len(re.findall(p, seq)) for p in FWD_PATS)
    rev = sum(len(re.findall(p, seq)) for p in REV_PATS)
    total = fwd + rev
    if total == 0:
        return 0.5, "unknown"
    frac = fwd / total
    if frac > 0.7:
        return frac, "coding-strand"
    elif frac < 0.3:
        return frac, "template-strand"
    else:
        return frac, "converging"


def find_longest_telo_run(seq, hexset=TELO_HEX):
    """Find the longest contiguous telomeric segment in seq."""
    n = len(seq)
    if n < 6:
        return 0, 0, 0
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

    # Find longest run of covered bases
    best_start = best_end = best_len = 0
    run_start = None
    for i in range(n):
        if covered[i]:
            if run_start is None:
                run_start = i
        else:
            if run_start is not None:
                run_len = i - run_start
                if run_len > best_len:
                    best_start = run_start
                    best_end = i
                    best_len = run_len
                run_start = None
    if run_start is not None:
        run_len = n - run_start
        if run_len > best_len:
            best_start = run_start
            best_end = n
            best_len = run_len

    return best_start, best_end, best_len


# ---------------------------------------------------------------------------
# Mode 1: VCF scanning
# ---------------------------------------------------------------------------
def scan_vcf(vcf_path, min_seq_len=30, min_telo_bp=18):
    """Scan VCF for insertions with telomeric content."""
    print(f"Scanning VCF: {vcf_path}", flush=True)

    results = []
    threshold_counts = {t: 0 for t in [0.10, 0.30, 0.50, 0.70, 0.85, 0.95]}
    n_variants = 0
    n_insertions = 0
    n_with_seq = 0

    opener = gzip.open if str(vcf_path).endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            n_variants += 1

            if n_variants % 100000 == 0:
                print(f"  Processed {n_variants:,} variants, "
                      f"{n_insertions:,} insertions, "
                      f"{len(results)} telo hits...", flush=True)

            fields = line.strip().split('\t', 8)  # only parse first 8 fields
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            sv_id = fields[2]
            ref = fields[3]
            alt = fields[4]

            # We want insertions — ALT longer than REF, or symbolic INS
            if alt.startswith('<'):
                # Symbolic allele — check INFO for SVTYPE=INS and SEQ
                if 'SVTYPE=INS' in fields[7] if len(fields) > 7 else '':
                    n_insertions += 1
                    # Try to extract sequence from INFO SEQ= field
                    info = fields[7]
                    seq_match = re.search(r'SEQ=([ACGTNacgtn]+)', info)
                    if seq_match:
                        ins_seq = seq_match.group(1).upper()
                    else:
                        continue
                else:
                    continue
            elif len(alt) > len(ref) + min_seq_len:
                # Explicit ALT sequence
                n_insertions += 1
                ins_seq = alt[len(ref):].upper()  # inserted portion
            elif ',' in alt:
                # Multi-allelic — check each
                for a in alt.split(','):
                    if len(a) > len(ref) + min_seq_len and not a.startswith('<'):
                        n_insertions += 1
                        ins_seq = a[len(ref):].upper()
                        break
                else:
                    continue
            else:
                continue

            n_with_seq += 1

            if len(ins_seq) < min_seq_len:
                continue

            # Check telomeric content
            cov = telo_coverage(ins_seq)

            for t in threshold_counts:
                if cov >= t:
                    threshold_counts[t] += 1

            if cov < 0.10:
                continue

            # Find longest telomeric run
            run_start, run_end, run_len = find_longest_telo_run(ins_seq)
            frac, orient = strand_ratio(ins_seq)

            result = {
                "chrom": chrom,
                "pos": pos,
                "sv_id": sv_id,
                "ins_length": len(ins_seq),
                "telo_coverage": round(cov, 4),
                "longest_telo_run": run_len,
                "run_start": run_start,
                "run_end": run_end,
                "orientation": orient,
                "fwd_frac": round(frac, 3),
                "seq_preview": ins_seq[:100],
            }

            # For high-purity hits, include more detail
            if cov >= 0.50:
                result["full_seq"] = ins_seq
                # Check for splice-like signals
                if len(ins_seq) >= 4:
                    result["first_4bp"] = ins_seq[:4]
                    result["last_4bp"] = ins_seq[-4:]

            results.append(result)

    print(f"\nVCF scan complete:", flush=True)
    print(f"  Total variants:        {n_variants:,}", flush=True)
    print(f"  Insertions:            {n_insertions:,}", flush=True)
    print(f"  With sequence ≥{min_seq_len}bp:  {n_with_seq:,}", flush=True)
    print(f"\nTelomeric insertions by threshold:", flush=True)
    for t in sorted(threshold_counts):
        print(f"  ≥{t*100:.0f}%: {threshold_counts[t]:,}", flush=True)

    return results, threshold_counts


# ---------------------------------------------------------------------------
# Mode 2: SVTig FASTA scanning
# ---------------------------------------------------------------------------
def scan_svtigs(fasta_path, min_len=30):
    """Scan SVTig FASTA for contigs with telomeric content."""
    print(f"Scanning SVTigs: {fasta_path}", flush=True)

    results = []
    threshold_counts = {t: 0 for t in [0.10, 0.30, 0.50, 0.70, 0.85, 0.95]}
    n_contigs = 0

    opener = gzip.open if str(fasta_path).endswith('.gz') else open

    curr_id = None
    parts = []

    def process_contig(cid, seq):
        nonlocal n_contigs
        n_contigs += 1
        if len(seq) < min_len:
            return

        cov = telo_coverage(seq)
        for t in threshold_counts:
            if cov >= t:
                threshold_counts[t] += 1

        if cov < 0.10:
            return

        run_start, run_end, run_len = find_longest_telo_run(seq)
        frac, orient = strand_ratio(seq)

        result = {
            "contig_id": cid,
            "length": len(seq),
            "telo_coverage": round(cov, 4),
            "longest_telo_run": run_len,
            "orientation": orient,
            "fwd_frac": round(frac, 3),
            "seq_preview": seq[:100],
        }
        if cov >= 0.50:
            result["full_seq"] = seq
        results.append(result)

    with opener(fasta_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if curr_id is not None:
                    process_contig(curr_id, ''.join(parts).upper())
                curr_id = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if curr_id is not None:
        process_contig(curr_id, ''.join(parts).upper())

    print(f"\nSVTig scan complete:", flush=True)
    print(f"  Total contigs:  {n_contigs:,}", flush=True)
    print(f"\nTelomeric contigs by threshold:", flush=True)
    for t in sorted(threshold_counts):
        print(f"  ≥{t*100:.0f}%: {threshold_counts[t]:,}", flush=True)

    return results, threshold_counts


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def write_report(results, threshold_counts, mode, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sort by coverage descending
    results.sort(key=lambda x: -x["telo_coverage"])

    # JSON details
    json_path = output_dir / f"telo_sv_hits_{mode}.json"
    with open(json_path, 'w') as f:
        json.dump({
            "mode": mode,
            "timestamp": datetime.now().isoformat(),
            "threshold_counts": {str(k): v for k, v in threshold_counts.items()},
            "n_hits_10pct": len(results),
            "hits": results[:500],  # cap at 500 for sanity
        }, f, indent=2)

    # TSV
    tsv_path = output_dir / f"telo_sv_hits_{mode}.tsv"
    with open(tsv_path, 'w') as f:
        if mode == "vcf":
            f.write("chrom\tpos\tsv_id\tins_length\ttelo_coverage\t"
                    "longest_telo_run\torientation\tfwd_frac\tseq_preview\n")
            for r in results:
                f.write(f"{r['chrom']}\t{r['pos']}\t{r.get('sv_id','')}\t"
                        f"{r['ins_length']}\t{r['telo_coverage']}\t"
                        f"{r['longest_telo_run']}\t{r['orientation']}\t"
                        f"{r['fwd_frac']}\t{r['seq_preview']}\n")
        else:
            f.write("contig_id\tlength\ttelo_coverage\tlongest_telo_run\t"
                    "orientation\tfwd_frac\tseq_preview\n")
            for r in results:
                f.write(f"{r['contig_id']}\t{r['length']}\t{r['telo_coverage']}\t"
                        f"{r['longest_telo_run']}\t{r['orientation']}\t"
                        f"{r['fwd_frac']}\t{r['seq_preview']}\n")

    # Markdown report
    md_path = output_dir / f"telo_sv_report_{mode}.md"
    with open(md_path, 'w') as f:
        f.write(f"# Telomeric SV Insertions — {mode.upper()} scan\n\n")
        f.write(f"**Date:** {datetime.now().isoformat()[:10]}\n\n")
        f.write("## Threshold counts\n\n")
        f.write("| Threshold | Count |\n|---|---:|\n")
        for t in sorted(threshold_counts):
            f.write(f"| ≥{t*100:.0f}% | {threshold_counts[t]:,} |\n")

        f.write(f"\n## Top hits (≥50% telomeric coverage)\n\n")
        hp = [r for r in results if r["telo_coverage"] >= 0.50]
        if not hp:
            f.write("None found.\n\n")
        else:
            for r in hp[:50]:
                if mode == "vcf":
                    f.write(f"### {r['chrom']}:{r['pos']} ({r['sv_id']})\n")
                    f.write(f"- Insertion length: {r['ins_length']} bp\n")
                else:
                    f.write(f"### {r['contig_id']}\n")
                    f.write(f"- Contig length: {r['length']} bp\n")
                f.write(f"- Telomeric coverage: {r['telo_coverage']:.1%}\n")
                f.write(f"- Longest telo run: {r['longest_telo_run']} bp\n")
                f.write(f"- Orientation: {r['orientation']} (fwd={r['fwd_frac']:.2f})\n")
                full = r.get('full_seq', r['seq_preview'])
                f.write(f"- Sequence: `{full[:200]}{'...' if len(full)>200 else ''}`\n\n")

    print(f"\nOutputs:", flush=True)
    print(f"  {json_path}", flush=True)
    print(f"  {tsv_path}", flush=True)
    print(f"  {md_path}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    output_dir = Path(__file__).parent / "results"

    if len(sys.argv) < 2:
        print("Usage:")
        print("  python scan_sv_telotrons.py vcf <path.vcf.gz>")
        print("  python scan_sv_telotrons.py svtig <path.fa.gz> [more.fa.gz ...]")
        sys.exit(1)

    mode = sys.argv[1]

    if mode == "vcf":
        vcf_path = sys.argv[2] if len(sys.argv) > 2 else "final-vcf.unphased.vcf.gz"
        results, counts = scan_vcf(vcf_path)
        write_report(results, counts, "vcf", output_dir)

    elif mode == "svtig":
        all_results = []
        all_counts = {t: 0 for t in [0.10, 0.30, 0.50, 0.70, 0.85, 0.95]}
        for fasta in sys.argv[2:]:
            results, counts = scan_svtigs(fasta)
            all_results.extend(results)
            for t in counts:
                all_counts[t] += counts[t]
        write_report(all_results, all_counts, "svtig", output_dir)

    else:
        print(f"Unknown mode: {mode}")
        sys.exit(1)

    # Print top hits
    print("\n" + "=" * 70)
    print("TOP TELOMERIC SV INSERTIONS")
    print("=" * 70)

    if mode == "vcf":
        top = sorted([r for r in (all_results if mode == "svtig" else results)],
                      key=lambda x: -x["telo_coverage"])
    else:
        top = sorted(all_results, key=lambda x: -x["telo_coverage"])

    for r in top[:20]:
        cov = r["telo_coverage"]
        if mode == "vcf":
            loc = f"{r['chrom']}:{r['pos']}"
            size = r["ins_length"]
        else:
            loc = r["contig_id"][:40]
            size = r["length"]
        print(f"  {cov:6.1%}  {size:>6}bp  {r['orientation']:>15s}  "
              f"telo_run={r['longest_telo_run']:>4}bp  {loc}")


if __name__ == "__main__":
    main()
