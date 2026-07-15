#!/usr/bin/env python
"""Linker segmentation of telotron / non-telo focal introns.

Parses every locus_text file under
  work/results/telotron_orthologs_v2/locus_text/{within_eimeria,outgroup}/{telotron_present,intron_present_nontelo}/*.txt

For each FOCAL intron sequence (the "-- intron --" line whose ortholog name ends with '*'):
 1. Strip lowercase / case markers, work in uppercase
 2. Find runs of >=3 tandem copies of any rotation / RC of {TTTAGGG, TTAGGG}
    allowing up to 1 mismatch per unit AND mean per-unit identity >= 0.85
 3. Everything between adjacent arrays (and intron-start <-> first / last <-> intron-end)
    is a LINKER
 4. Architecture = ordered string like L|A+|L|A-|L
 5. Emit JSONL (per locus) + per-linker TSV + per-locus TSV (one row per UNIQUE
    locus_text filename; panel and category are emitted as semicolon-joined lists)
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import SPECIES_BY_ACC  # single source of truth for the focal panel

# ROOT / OUT_DIR are resolved at runtime from CLI args / env vars / repo-relative
# defaults so the script no longer hardcodes absolute /scratch1 paths (commit
# 21fd9d4 "Consolidate pipeline data dirs under work/"; FIXES_PRIORITIZED.md #103).
_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ROOT = Path(
    os.environ.get(
        "LINKER_SEGMENTATION_ROOT",
        _REPO_ROOT / "work" / "results" / "telotron_orthologs_v2" / "locus_text",
    )
)
DEFAULT_OUT_DIR = Path(
    os.environ.get(
        "LINKER_SEGMENTATION_OUT_DIR",
        _REPO_ROOT / "work" / "results" / "mechanism_deepdive",
    )
)
# Module-level aliases preserved for backwards compatibility with any importer
# that does `from linker_segmentation import ROOT, OUT_DIR`. main() resolves the
# effective values from argparse and rebinds these globals before use.
ROOT = DEFAULT_ROOT
OUT_DIR = DEFAULT_OUT_DIR

PANELS = ("within_eimeria", "outgroup")
CATEGORIES = ("telotron_present", "intron_present_nontelo")

# Re-export under the legacy local name so existing consumers
# (e.g. mask_telotron_arrays imports GCF_TO_SPECIES from here) keep working.
GCF_TO_SPECIES = SPECIES_BY_ACC

# --- motif set + rotations + RCs ---
# Restricted to the Apicomplexan canonical (TTTAGGG, per config.yaml:canonical_telomere_motifs)
# plus the vertebrate fallback (TTAGGG). The 5-bp arthropod motif TTAGG and the
# 6-bp ciliate motif TTGGGG were dropped: they false-positively segment AT-rich
# linker noise at min_copies=3 + 1 mm/unit, inflating the n_arrays>=6 tail and
# splitting real linkers around fake arrays (verify_linker_architecture HIGH;
# telotron_mechanism_deepdive_2026-06-03 segmenter audit; FIXES_PRIORITIZED.md #9).
BASE_MOTIFS = ("TTTAGGG", "TTAGGG")

def rc(s: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp[b] for b in reversed(s))

def all_rotations(s: str):
    return {s[i:] + s[:i] for i in range(len(s))}

def build_motif_table():
    """Each entry: (label, length, orient, frozenset(units))
    where the units are all rotations (G-rich keeps original orient '+',
    RC of G-rich is the C-rich strand '-').
    """
    table = []
    seen = {}
    for m in BASE_MOTIFS:
        plus = all_rotations(m)
        minus = all_rotations(rc(m))
        # use canonical label = m for '+', rc(m) for '-'
        table.append((m, len(m), "+", frozenset(plus)))
        table.append((rc(m), len(m), "-", frozenset(minus)))
    return table

MOTIF_TABLE = build_motif_table()

# Identity floors per TRF (Benson 1999 NAR 27:573) / tidk (Brown 2025
# Bioinformatics btaf049) / Azzalin 2001 ITS criterion (<1 mm/unit ≡ ≥6/7
# identity). Defaults: <=1 mm per unit AND >=85% mean unit identity across the
# extended run. Verifier: review_linker_architecture HIGH find_arrays; new
# concern NC-1 aggregate-identity-floor. FIXES_PRIORITIZED.md #15.
MAX_MM_PER_UNIT = 1
MIN_MEAN_UNIT_IDENT = 0.85

def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)

def best_unit_mm(window: str, units: frozenset) -> int:
    """Minimum Hamming distance from `window` to any unit in `units`.
    Returns len(window)+1 if `units` is empty or no unit shares the window length.
    """
    if window in units:
        return 0
    best = len(window) + 1
    for u in units:
        d = hamming(window, u)
        if d < best:
            best = d
            if best == 0:
                return 0
    return best

def match_unit(window: str, units: frozenset, max_mm: int = MAX_MM_PER_UNIT) -> bool:
    return best_unit_mm(window, units) <= max_mm

def find_arrays(seq: str):
    """Greedy non-overlapping scan: for each position, try each motif/orient,
    extend while consecutive non-overlapping units match (<=MAX_MM_PER_UNIT);
    if >=3 copies AND mean per-unit identity >= MIN_MEAN_UNIT_IDENT, record
    the array and jump past it.

    Per-copy identity is enforced by `match_unit` (<=1 mm per unit); aggregate
    identity is enforced by the mean-mm/ulen check below. "Best motif" tiebreaks
    on (span, fewer total mm, longer unit) so a 5-bp slop-tolerant motif cannot
    silently outscore a 7-bp perfect-match motif (verify_linker_architecture
    new concern NC-1).

    Returns list of dicts: {start, end, motif, orient, copies, unit_len, mm_total, mean_ident}.
    """
    arrays = []
    n = len(seq)
    i = 0
    while i < n:
        best = None  # (end, copies, motif_label, orient, unit_len, mm_total)
        for label, ulen, orient, units in MOTIF_TABLE:
            if i + ulen > n:
                continue
            # must start with a match
            first_mm = best_unit_mm(seq[i:i+ulen], units)
            if first_mm > MAX_MM_PER_UNIT:
                continue
            j = i
            copies = 0
            mm_total = 0
            while j + ulen <= n:
                window_mm = best_unit_mm(seq[j:j+ulen], units)
                if window_mm > MAX_MM_PER_UNIT:
                    break
                copies += 1
                mm_total += window_mm
                j += ulen
            if copies < 3:
                continue
            # aggregate identity floor: mean per-unit identity over the run
            mean_ident = 1.0 - (mm_total / float(copies * ulen)) if copies > 0 else 0.0
            if mean_ident < MIN_MEAN_UNIT_IDENT:
                continue
            span = j - i
            # tiebreak: (longer span, fewer mm_total, longer unit_len)
            key = (span, -mm_total, ulen)
            if best is None or key > (best[0] - i, -best[5], best[4]):
                best = (j, copies, label, orient, ulen, mm_total)
        if best is not None:
            end, copies, label, orient, ulen, mm_total = best
            mean_ident = 1.0 - (mm_total / float(copies * ulen))
            arrays.append({
                "start": i,
                "end": end,
                "motif": label,
                "orient": orient,
                "copies": copies,
                "unit_len": ulen,
                "mm_total": mm_total,
                "mean_ident": round(mean_ident, 4),
            })
            i = end
        else:
            i += 1
    # NB: a previous "merge overlapping arrays" pass lived here but was dead
    # code — `i = end` in the scan loop prevents any overlap from arising
    # (verify_linker_architecture LOW; FIXES_PRIORITIZED.md #102). Removed.
    return arrays

def gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    g = sum(1 for b in seq if b in "GC")
    return round(100.0 * g / len(seq), 2)

# ---- parsing the locus_text files ----

COLOR_RE = re.compile(r"<color>\s*(.*?)\s*</color>", re.S)

def parse_focal_intron(path: Path) -> str | None:
    """Return uppercase intron DNA for the focal species (the line whose ortho
    name ends with '*' inside the UNALIGNED DNA `-- intron --` section).
    """
    try:
        text = path.read_text(errors="replace")
    except Exception:
        return None
    # find the UNALIGNED DNA section
    m = re.search(r"### UNALIGNED — DNA ###(.*?)(?:###|\Z)", text, re.S)
    if not m:
        return None
    section = m.group(1)
    # find the -- intron -- block (up to next "-- ")
    m2 = re.search(r"-- intron --\s*(.*?)(?:\n-- |\Z)", section, re.S)
    if not m2:
        return None
    block = m2.group(1)
    # within block, find the focal line: starts at line beginning, ortho id has '*'
    for line in block.splitlines():
        line_strip = line.strip()
        if not line_strip:
            continue
        # focal line begins with something like 'GCF_000499385.1*'
        # take any line where the first token ends with '*'
        first = line_strip.split()[0]
        if first.endswith("*"):
            cm = COLOR_RE.search(line)
            if cm:
                return cm.group(1).upper().replace(" ", "").replace("-", "")
    return None

def parse_focal_gcf(path: Path) -> str | None:
    """Parse the focal GCF from filename: GCF_xxxxx.1__contig_start-end.txt"""
    name = path.name
    if "__" in name:
        gcf = name.split("__", 1)[0]
        return gcf
    return None

def telotron_id_from_path(path: Path) -> str:
    return path.stem  # GCF_000499385.1__NW_013651667.1_11188-11409

def process_one(args):
    panel, category, path_str = args
    path = Path(path_str)
    tid = telotron_id_from_path(path)
    gcf = parse_focal_gcf(path)
    species = GCF_TO_SPECIES.get(gcf, gcf or "unknown")
    intron = parse_focal_intron(path)
    if intron is None or len(intron) == 0:
        return None
    # arrays + linkers
    arrays = find_arrays(intron)
    n = len(intron)
    # build linker list interleaved with arrays
    linkers = []
    cursor = 0
    for k, a in enumerate(arrays):
        if a["start"] > cursor:
            linker_seq = intron[cursor:a["start"]]
            if k == 0:
                pos = "5p"
            else:
                pos = "mid"
            linkers.append({
                "position": pos,
                "len": len(linker_seq),
                "gc": gc_pct(linker_seq),
                "seq": linker_seq,
                "prev_motif": arrays[k-1]["motif"] if k > 0 else None,
                "prev_orient": arrays[k-1]["orient"] if k > 0 else None,
                "next_motif": a["motif"],
                "next_orient": a["orient"],
                "idx": k,  # index of next array
            })
        cursor = a["end"]
    # trailing 3' linker (or whole-intron linker if no arrays)
    if cursor < n:
        if not arrays:
            pos = "whole"
            prev_motif = prev_orient = next_motif = next_orient = None
        else:
            pos = "3p"
            prev_motif = arrays[-1]["motif"]
            prev_orient = arrays[-1]["orient"]
            next_motif = next_orient = None
        linker_seq = intron[cursor:n]
        linkers.append({
            "position": pos,
            "len": len(linker_seq),
            "gc": gc_pct(linker_seq),
            "seq": linker_seq,
            "prev_motif": prev_motif,
            "prev_orient": prev_orient,
            "next_motif": next_motif,
            "next_orient": next_orient,
            "idx": len(arrays),
        })
    # also: case where there are arrays but cursor==n (no trailing linker -> intron ends in array, but 3' GT/AG: usually intron ends in NNAG so there's always a tail).
    # architecture string interleaving
    arch_parts = []
    for k, a in enumerate(arrays):
        # leading linker if k==0 and arrays[0].start>0
        if k == 0 and a["start"] > 0:
            arch_parts.append("L")
        elif k > 0:
            arch_parts.append("L")
        arch_parts.append(f"A{a['orient']}")
    if not arrays:
        arch_parts.append("L")
    else:
        # trailing linker
        if arrays[-1]["end"] < n:
            arch_parts.append("L")
    architecture = "|".join(arch_parts)
    total_array_len = sum(a["end"] - a["start"] for a in arrays)
    total_linker_len = sum(l["len"] for l in linkers)
    rec = {
        "telotron": tid,
        "species": species,
        "panel": panel,
        "category": category,
        "intron_len": n,
        "n_arrays": len(arrays),
        "architecture": architecture,
        "total_array_len": total_array_len,
        "total_linker_len": total_linker_len,
        "arrays": arrays,
        "linkers": linkers,
    }
    return rec

def main():
    ap = argparse.ArgumentParser(
        description="Segment focal-intron sequences into arrays + linkers."
    )
    ap.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="locus_text directory root (overrides $LINKER_SEGMENTATION_ROOT)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="output directory (overrides $LINKER_SEGMENTATION_OUT_DIR)",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("LINKER_SEGMENTATION_WORKERS", "32")),
    )
    args = ap.parse_args()

    # Rebind module-level globals so any consumer that imports ROOT / OUT_DIR
    # after main() runs sees the resolved values.
    global ROOT, OUT_DIR
    ROOT = args.root
    OUT_DIR = args.out_dir
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Collect UNIQUE locus_text files. The same telotron filename can live under
    # multiple (panel, category) subdirs; emitting one row per (panel, category)
    # × locus inflates architecture_per_locus.tsv ~1.65× (verifier: 742 rows for
    # 451 unique loci). Iterate via rglob and dedupe by filename; record the
    # full panel/category list as semicolon-joined columns
    # (FIXES_PRIORITIZED.md #4/#16; review_linker_architecture HIGH).
    by_name: dict[str, dict] = {}
    for p in sorted(ROOT.rglob("*.txt")):
        rel = p.relative_to(ROOT)
        parts = rel.parts
        panel = parts[0] if len(parts) >= 1 and parts[0] in PANELS else "unknown"
        cat = parts[1] if len(parts) >= 2 and parts[1] in CATEGORIES else "unknown"
        rec = by_name.setdefault(
            p.name,
            {"path": str(p), "panels": set(), "categories": set()},
        )
        rec["panels"].add(panel)
        rec["categories"].add(cat)
        # Keep the first path; all copies of a given filename should have
        # identical content (locus_text files are partitioned, not edited).
    tasks = []
    for name, rec in sorted(by_name.items()):
        panel_str = ";".join(sorted(rec["panels"]))
        cat_str = ";".join(sorted(rec["categories"]))
        tasks.append((panel_str, cat_str, rec["path"]))
    print(
        f"[linker_segmentation] {len(tasks)} unique locus files to process "
        f"(from {sum(1 for _ in ROOT.rglob('*.txt'))} total entries under {ROOT})",
        file=sys.stderr,
    )

    results = []
    n_skipped = 0
    skipped_paths: list[str] = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(process_one, t): t for t in tasks}
        for f in as_completed(futures):
            r = f.result()
            if r is None:
                n_skipped += 1
                skipped_paths.append(futures[f][2])
            else:
                results.append(r)
    print(
        f"[linker_segmentation] {len(results)} non-empty loci parsed; "
        f"{n_skipped} skipped (intron parse failed or empty)",
        file=sys.stderr,
    )
    if n_skipped and len(skipped_paths) <= 20:
        for sp in skipped_paths:
            print(f"[linker_segmentation]   skipped: {sp}", file=sys.stderr)

    # outputs
    jsonl_path = OUT_DIR / "linker_segmentation.jsonl"
    per_linker_path = OUT_DIR / "linker_segmentation.tsv"
    per_locus_path = OUT_DIR / "architecture_per_locus.tsv"

    with jsonl_path.open("w") as jf, per_linker_path.open("w") as lf, per_locus_path.open("w") as af:
        lf.write("telotron\tspecies\tpanel\tcategory\tposition\tlen\tgc\tseq\tprev_motif\tprev_orient\tnext_motif\tnext_orient\n")
        af.write("telotron\tspecies\tpanel\tcategory\tn_arrays\tarchitecture\ttotal_array_len\ttotal_linker_len\tintron_len\n")
        for r in results:
            jf.write(json.dumps(r) + "\n")
            af.write("\t".join([
                r["telotron"], r["species"], r["panel"], r["category"],
                str(r["n_arrays"]), r["architecture"],
                str(r["total_array_len"]), str(r["total_linker_len"]),
                str(r["intron_len"]),
            ]) + "\n")
            for l in r["linkers"]:
                lf.write("\t".join([
                    r["telotron"], r["species"], r["panel"], r["category"],
                    l["position"], str(l["len"]), str(l["gc"]), l["seq"],
                    str(l.get("prev_motif") or ""),
                    str(l.get("prev_orient") or ""),
                    str(l.get("next_motif") or ""),
                    str(l.get("next_orient") or ""),
                ]) + "\n")

    # summary stats
    import statistics as st
    from collections import Counter
    n_loci = len(results)
    all_linkers = [l for r in results for l in r["linkers"]]
    n_linkers = len(all_linkers)
    lengths = [l["len"] for l in all_linkers]
    gcs = [l["gc"] for l in all_linkers if l["len"] > 0]
    def pct(arr, q):
        if not arr:
            return None
        arr_sorted = sorted(arr)
        k = max(0, min(len(arr_sorted)-1, int(round((len(arr_sorted)-1) * q / 100.0))))
        return arr_sorted[k]
    arch_counter = Counter(r["architecture"] for r in results)
    top_arch = arch_counter.most_common(20)
    summary = {
        "n_loci_processed": n_loci,
        "n_total_linker_segments": n_linkers,
        "linker_len_mean": round(st.mean(lengths), 2) if lengths else None,
        "linker_len_median": st.median(lengths) if lengths else None,
        "linker_len_p90": pct(lengths, 90),
        "linker_gc_mean": round(st.mean(gcs), 2) if gcs else None,
        "linker_gc_p10": pct(gcs, 10),
        "linker_gc_p90": pct(gcs, 90),
        "top_architectures": top_arch,
    }
    print(json.dumps(summary, indent=2, default=str))

if __name__ == "__main__":
    main()
