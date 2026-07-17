#!/usr/bin/env python3
"""Mask every telotron intron into a G/A/L architecture cartoon, align the
cartoons on the inverted-repeat junction, and present the result locus_text-style.

Per focal telotron intron (work/results/telotron_orthologs_v2/locus_text/*/telotron_present/*.txt):
  - find telomeric-repeat arrays with linker_segmentation.find_arrays
  - mask each base:  G = G-rich array (orient '+', TTTAGGG side)
                     A = A-rich array (orient '-', CCCTAAA side)
                     L = linker (everything not in an array)

Alignment: telotrons are bidirectional (A-rich arm … linker … G-rich arm). MAFFT
on a 3-letter low-complexity cartoon just scatters gaps, so instead we orient every
cartoon A→G (reverse-complementing the minority that run G→A — the arms are mutual
revcomps) and anchor on the first inter-array A→G transition coordinate, padding
left/right. Loci with only one polarity (no A-arm or no G-arm) have no biological
A↔G meeting point and are flagged non_anchorable / left-padded.

CAVEAT: the G/A/L cartoon is low-complexity by construction; per-column consensus
statistics on this MSA are NOT base-orthology statements (Frith 2011 PMC3327517;
Bzikadze & Pevzner 2022). Use the MSA for structural visualisation only.

Outputs (work/results/telotron_masked/):
  telotron_masked.fasta       one G/A/L string per telotron (gene-sense)
  telotron_masked.msa.fasta   junction-anchored, A→G-oriented alignment (gaps = '-')
  telotron_masked_msa.txt     locus_text-style view (rows wrapped in <color>, ruler)
  telotron_masked_msa.html    self-contained colored view (A green · G gold · L slate)
"""
from __future__ import annotations
import html as _html
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import SPECIES_ABBREV as SP_ABBREV, SPECIES_BY_ACC
from linker_segmentation import (  # reuse the audited segmentation
    find_arrays, parse_focal_intron, parse_focal_gcf, telotron_id_from_path,
    ROOT, all_rotations,
)
# Single source of truth for the focal-panel GCF -> species map is
# _common.SPECIES_BY_ACC (config.yaml-derived); do not import the legacy
# GCF_TO_SPECIES alias from linker_segmentation.
GCF_TO_SPECIES = SPECIES_BY_ACC

PANELS = ("within_eimeria", "outgroup")
CATEGORY = "telotron_present"
OUT_DIR = Path("work/results/telotron_masked")
CLEAN_CHAR = {"+": "G", "-": "A"}           # perfect G-rich / A-rich repeat unit
DEGEN_CHAR = {"+": "Ğ", "-": "Ǎ"}           # degenerate (<=1 mismatch) repeat unit
GSET, ASET = set("GĞ"), set("AǍ")           # arm membership (any purity)
SWAP = {"A": "G", "G": "A", "Ǎ": "Ğ", "Ğ": "Ǎ", "L": "L"}  # arm flip (revcomp)
LINKER_CHAR = "L"


def mask_intron(seq: str) -> tuple[str, str, list]:
    """Return (masked string, architecture, arrays). Mask alphabet per base:
       A/G = perfect A-rich/G-rich repeat unit, Ǎ/Ğ = degenerate (<=1 mm) unit,
       L = linker (anything not inside a detected >=3-copy array)."""
    arrays = find_arrays(seq)
    mask = [LINKER_CHAR] * len(seq)
    for a in arrays:
        ulen = a["unit_len"]
        perfect = all_rotations(a["motif"])         # exact rotations of this motif
        clean, degen = CLEAN_CHAR[a["orient"]], DEGEN_CHAR[a["orient"]]
        p = a["start"]
        while p < a["end"]:
            ch = clean if seq[p:p + ulen] in perfect else degen
            for q in range(p, min(p + ulen, a["end"])):
                mask[q] = ch
            p += ulen
    parts, n = [], len(seq)
    if arrays:
        if arrays[0]["start"] > 0:
            parts.append("L")
        for k, a in enumerate(arrays):
            if k > 0:
                parts.append("L")
            parts.append("A" + a["orient"])
        if arrays[-1]["end"] < n:
            parts.append("L")
    else:
        parts.append("L")
    return "".join(mask), "|".join(parts), arrays


def canonicalize(mask: str) -> tuple[str, str]:
    """Return (mask unchanged, orientation) where orientation is 'A_first',
    'G_first', or 'single'.

    THIS NO LONGER FLIPS, because flipping cannot work — arm order is
    revcomp-INVARIANT. The old version reverse-complemented (reverse + swap
    A/G, since the arms are mutual revcomps) whenever G mass sat 5' of A mass,
    intending to normalise every locus to A-rich-5'/G-rich-3'. That is
    algebraically impossible: revcomp([G-arm][linker][A-arm]) is
    [revcomp(A-arm)][revcomp(linker)][revcomp(G-arm)] = [G-arm][linker][A-arm]
    again, because revcomp of an A-arm IS a G-arm. Formally, the flip sends
    mean(gi) -> n-1-mean(ai) and mean(ai) -> n-1-mean(gi), so the predicate
    mean(gi) < mean(ai) that triggered the flip still holds afterwards.
    Verified: a 'GGGGGGGGGGGGGGLLAAAAAAAAAAAAAA' cartoon came back
    byte-identical while reporting flipped=True.

    The consequence was severe: structural_junction() only looked for an A->G
    transition, so every G-first bidirectional locus (GT-F-R-AG is ~44% of
    loci) returned None and was mislabeled non_anchorable=single_polarity —
    while its own polarities field reported both arms present. Those loci were
    dropped from the junction-anchored MSA and inflated the reported
    single-polarity count.

    Orientation is now REPORTED rather than normalised away, which is also the
    scientifically right call: dyad polarity is a result in its own right
    (~99% of Mixed loci are G-rich-5'), not a strand-convention artifact to be
    canonicalised out.
    """
    gi = [i for i, c in enumerate(mask) if c in GSET]
    ai = [i for i, c in enumerate(mask) if c in ASET]
    if not (gi and ai):
        return mask, "single"
    return mask, ("G_first" if sum(gi) / len(gi) < sum(ai) / len(ai) else "A_first")


def structural_junction(mask: str) -> int | None:
    """Coordinate of the first inter-array polarity transition, in EITHER
    direction (A-arm -> G-arm or G-arm -> A-arm).

    Orientation-agnostic by necessity: canonicalize() cannot reorient a
    G-first locus (see above), so a junction finder that only recognises A->G
    silently discards the entire G-first class. Returns None only for
    genuinely single-polarity loci.
    """
    has_a = any(c in ASET for c in mask)
    has_g = any(c in GSET for c in mask)
    if not (has_a and has_g):
        return None
    last_arm = None  # 'A' or 'G'
    for i, c in enumerate(mask):
        if c in ASET:
            if last_arm == "G":
                return i           # first G->A transition
            last_arm = "A"
        elif c in GSET:
            if last_arm == "A":
                return i           # first A->G transition
            last_arm = "G"
    # Both polarities present but no transition found is unreachable: the first
    # arm character sets last_arm, and the first character of the other
    # polarity necessarily triggers a return above.
    raise AssertionError(
        f"both polarities present but no transition found in cartoon: {mask!r}"
    )


def collect():
    recs = []
    for panel in PANELS:
        d = ROOT / panel / CATEGORY
        if not d.is_dir():
            continue
        for p in sorted(d.glob("*.txt")):
            intron = parse_focal_intron(p)
            if not intron:
                continue
            masked, arch, arrays = mask_intron(intron)
            canon, orientation = canonicalize(masked)   # canon == masked; no flip is possible
            jx = structural_junction(canon)             # None ONLY for single-polarity loci
            orients = {a["orient"] for a in arrays}
            gcf = parse_focal_gcf(p)
            recs.append({
                "id": telotron_id_from_path(p),
                "species": GCF_TO_SPECIES.get(gcf, gcf or "unknown"),
                "len": len(intron),
                "masked": masked,      # gene-sense
                "canon": canon,        # gene-sense (orientation reported, not normalised)
                "arch": arch,
                "orientation": orientation,     # 'A_first' | 'G_first' | 'single'
                "anchorable": jx is not None,   # False = genuinely single-polarity
                "junction": jx,                 # first polarity transition, either direction
                "n_arrays": len(arrays),
                "polarities": "".join(sorted(orients)) or "none",
            })
    return recs


def anchor_align(recs):
    """Junction-anchored alignment of the A→G-oriented cartoons.

    Junction = first inter-array A→G transition coordinate (`r["junction"]`,
    set by ``structural_junction``). Single-polarity loci (``anchorable=False``)
    have no biological A↔G meeting point; they are left-justified at the global
    junction column with `r["junction"]` set to ``len(canon)`` so the row pads
    only to the right (the prior behaviour of "first G" anchored these loci at
    column 0 — a degenerate anchor, replaced here)."""
    for r in recs:
        if r["anchorable"]:
            j = r["junction"]
        else:
            # Non-anchorable: place the whole cartoon left of the junction
            # column so anchored loci still align cleanly.
            j = len(r["canon"])
        r["left"], r["right"] = r["canon"][:j], r["canon"][j:]
    maxL = max(len(r["left"]) for r in recs)
    maxR = max((len(r["right"]) for r in recs), default=0)
    for r in recs:
        r["aln"] = r["left"].rjust(maxL, "-") + r["right"].ljust(maxR, "-")
    return maxL, maxL + maxR


def write_fasta(recs, path, key, note):
    with path.open("w") as fh:
        for r in recs:
            anchor_tag = " non_anchorable=single_polarity" if not r["anchorable"] else ""
            fh.write(f">{r['id']} species={r['species']} len={r['len']} "
                     f"arch={r['arch']} polarities={r['polarities']} "
                     f"orientation={r['orientation']}{anchor_tag}{note}\n")
            fh.write(r[key] + "\n")


def label(r):
    sp = SP_ABBREV.get(r["species"], r["species"][:4])
    coords = r["id"].split("__", 1)[1] if "__" in r["id"] else r["id"]
    return f"{sp}:{coords}"


def order(recs):
    return sorted(recs, key=lambda r: (r["arch"].count("A"), r["arch"], len(r["left"]), r["species"]))


def present_txt(recs, junction, width, out_txt):
    rows = order(recs)
    n_nonanchor = sum(1 for r in rows if not r["anchorable"])
    labels = [label(r) for r in rows]
    idw = min(46, max(len(x) for x in labels))
    OPEN, CLOSE = "<color>", "</color>"
    pad = idw + 1 + len(OPEN) + 1
    tens = "".join(str((c // 10) % 10) if (c + 1) % 10 == 0 else " " for c in range(width))
    L = [
        "=" * 78,
        f"TELOTRON ARRAY-ARCHITECTURE MSA   (n={len(rows)} telotrons; "
        f"{n_nonanchor} single-polarity / non-anchorable, left-padded)",
        "mask:  A/Ǎ = A-rich array perfect/degenerate (CCCTAAA side) · "
        "G/Ğ = G-rich perfect/degenerate (TTTAGGG side) · L = linker (non-repeat)",
        "oriented A-rich(5') -> G-rich(3'); anchored on the first inter-array",
        "A->G transition coordinate (caret below). Loci with only one polarity",
        "are flagged non_anchorable and contribute padding, not a structural call.",
        "CAVEAT: the G/A/L cartoon is low-complexity by construction; per-column",
        "consensus / column-wise statistics on this MSA are NOT base-orthology",
        "statements (cf. Frith 2011 PMC3327517; Bzikadze & Pevzner 2022).",
        "=" * 78, "",
        "### ALIGNED — array cartoon ###", "",
        " " * pad + tens,
        " " * (pad + junction) + "^  junction (A-arm | G-arm)",
    ]
    for r, lab in zip(rows, labels):
        tag = "  [non-anchorable]" if not r["anchorable"] else ""
        L.append(f"{lab:<{idw}.{idw}} {OPEN} {r['aln']} {CLOSE}{tag}")
    out_txt.write_text("\n".join(L) + "\n")


HTML_TMPL = """<!DOCTYPE html><html><head><meta charset="utf-8"><title>Telotron array MSA</title>
<style>
 body{{background:#1e2127;color:#c9d1d9;font-family:-apple-system,sans-serif;margin:0}}
 header{{padding:10px 16px;border-bottom:1px solid #363b45;position:sticky;top:0;background:#23262d;z-index:2}}
 h1{{font-size:14px;margin:0 0 4px}} .key{{font-size:12px;color:#8b949e}}
 .caveat{{font-size:11px;color:#c98b4a;margin-top:4px}}
 .A{{color:#4eb87a;font-weight:600}} .a2{{color:#2f6b48}}
 .G{{color:#e0b24a;font-weight:600}} .g2{{color:#8a6f33}}
 .L{{color:#7b8694}} .gap{{color:#2b2f38}}
 .nax{{color:#c98b4a;font-weight:600}}
 .key b.A{{color:#4eb87a}} .key b.G{{color:#e0b24a}} .key b.L{{color:#7b8694}}
 .wrap{{overflow:auto;padding:10px 16px}}
 pre{{margin:0;font:12px/1.35 ui-monospace,Menlo,Consolas,monospace;white-space:pre}}
 .row{{white-space:pre}} .lab{{color:#8b949e}} .jx{{border-left:1px solid #e0625b55}}
</style></head><body>
<header><h1>Telotron array-architecture MSA (n={n}; {n_nonanchor} single-polarity / non-anchorable)</h1>
<div class="key">mask: <b class="A">A</b>/<span class="a2">Ǎ</span> A-rich array perfect/degenerate (CCCTAAA) ·
 <b class="G">G</b>/<span class="g2">Ğ</span> G-rich perfect/degenerate (TTTAGGG) ·
 <b class="L">L</b> linker (non-repeat) —
 oriented A-rich(5')→G-rich(3'), anchored on the first inter-array A→G transition (col {jx});
 single-polarity loci are flagged <span class="nax">[non-anchorable]</span> and left-padded.</div>
<div class="caveat">CAVEAT: the G/A/L cartoon is low-complexity by construction; per-column
 consensus / column-wise statistics on this MSA are NOT base-orthology statements
 (Frith 2011 PMC3327517; Bzikadze &amp; Pevzner 2022).</div></header>
<div class="wrap"><pre>{body}</pre></div></body></html>"""


def present_html(recs, junction, out_html):
    rows = order(recs)
    n_nonanchor = sum(1 for r in rows if not r["anchorable"])
    idw = min(46, max(len(label(r)) for r in rows))
    cls = {"A": "A", "G": "G", "L": "L", "-": "gap", "Ǎ": "a2", "Ğ": "g2"}
    body = []
    for r in rows:
        lab = _html.escape(f"{label(r):<{idw}.{idw}}")
        spans, run, rc = [], [], None
        for c in r["aln"]:
            k = cls[c]
            if k != rc and run:
                spans.append(f'<span class="{rc}">{"".join(run)}</span>'); run = []
            rc = k; run.append(c)
        if run:
            spans.append(f'<span class="{rc}">{"".join(run)}</span>')
        tag = '  <span class="nax">[non-anchorable]</span>' if not r["anchorable"] else ""
        body.append(f'<span class="lab">{lab}</span>  ' + "".join(spans) + tag)
    out_html.write_text(HTML_TMPL.format(n=len(rows), n_nonanchor=n_nonanchor,
                                         jx=junction, body="\n".join(body)))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    recs = collect()
    if not recs:
        sys.exit("no telotron introns found")
    write_fasta(recs, OUT_DIR / "telotron_masked.fasta", "masked", "")
    junction, width = anchor_align(recs)
    write_fasta(recs, OUT_DIR / "telotron_masked.msa.fasta", "aln", " anchored")
    present_txt(recs, junction, width, OUT_DIR / "telotron_masked_msa.txt")
    present_html(recs, junction, OUT_DIR / "telotron_masked_msa.html")

    n_nonanchor = sum(1 for r in recs if not r["anchorable"])
    print(f"{len(recs)} telotrons | alignment width {width} (junction at col {junction}) "
          f"| {n_nonanchor} non-anchorable (single-polarity)")
    print("top architectures:")
    for a, c in Counter(r["arch"] for r in recs).most_common(12):
        print(f"   {c:4d}  {a}")
    # Dyad orientation is reported, not normalised: no revcomp can change arm
    # order (see canonicalize). The old "flipped to canonical orientation: N"
    # counter was meaningless — the flip returned the input unchanged.
    print("dyad orientation:", dict(Counter(r["orientation"] for r in recs)))
    for f in ("telotron_masked.fasta", "telotron_masked.msa.fasta",
              "telotron_masked_msa.txt", "telotron_masked_msa.html"):
        print("wrote", OUT_DIR / f)


if __name__ == "__main__":
    main()
