#!/usr/bin/env python3
"""Per-telotron-locus figure (compiled multipage PDF).

For each telotron, one page with two panels:
  TOP  protein MSA of the flanking exons, zoomed to +/- WINDOW residues around the
       intron insertion point (focal gene + orthologs); cells coloured by match to the
       focal residue, gaps white, premature stops (frameshift/pseudogene) red. The
       intron position is marked between the left- and right-flank exon residues.
  BOTTOM telotron-vs-orthologous-intron DNA alignment as telomere-highlighted tracks
       (telomeric-repeat bases orange, other bases grey, gaps white), one row per
       sequence, annotated with telomere-fraction and any FRAMESHIFT@junction flag.

Reads the outputs of telotron_ortholog_align.py (protein_msa/, intron_dna/,
ortholog_map.tsv). Writes <outdir>/telotron_ortholog_loci.pdf (+ optional PNGs).
"""
import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import rc  # noqa: E402


def read_aln(path):
    """Return list of (name, seq) preserving order."""
    out, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def telomeric_mask(seq, motif):
    """Boolean per (non-gap) position: is this base part of a >=2-unit telomeric run?
    Aligns to the best phase of motif or its revcomp over each maximal ungapped block."""
    units = [motif, rc(motif)]
    L = len(motif)
    mask = [False] * len(seq)
    # operate on ungapped index mapping
    ung = [i for i, c in enumerate(seq) if c != "-"]
    s = "".join(seq[i] for i in ung).upper()
    if len(s) < L:
        return mask
    best = None
    bestscore = -1
    for u in units:
        for ph in range(L):
            patt = (u[ph:] + u * (len(s) // L + 2))[:len(s)]
            sc = sum(1 for a, b in zip(s, patt) if a == b)
            if sc > bestscore:
                bestscore, best = sc, patt
    # mark positions matching the best telomeric phase, but only inside runs of >=10 matched
    matched = [a == b for a, b in zip(s, best)]
    i = 0
    while i < len(matched):
        if matched[i]:
            j = i
            while j < len(matched) and (matched[j] or (j - i < 14)):  # bridge <=14bp gaps in run
                j += 1
            if sum(matched[i:j]) >= 8:
                for k in range(i, j):
                    mask[ung[k]] = matched[k]
            i = j
        else:
            i += 1
    return mask


def focal_col(focal_seq, residue):
    """Map an ungapped focal residue index to its MSA column index."""
    n = 0
    for col, c in enumerate(focal_seq):
        if c != "-":
            n += 1
            if n == residue:
                return col
    return len(focal_seq) - 1


def short(name):
    """Trim alignment header to a compact row label."""
    gid = name.split("|")[0]
    tags = []
    for m in re.finditer(r"(id|fs|stop|telo)=([0-9.]+)", name):
        if m.group(1) == "id" and float(m.group(2)) > 0:
            tags.append(f"id{m.group(2)}")
        elif m.group(1) in ("fs", "stop") and float(m.group(2)) > 0:
            tags.append(f"{m.group(1)}{int(float(m.group(2)))}")
        elif m.group(1) == "telo":
            tags.append(f"telo{m.group(2)}")
    if "FOCAL" in name:
        tags.insert(0, "FOCAL")
    if "FRAMESHIFT@junction" in name:
        tags.append("FS@junc")
    return gid + (" " + ",".join(tags) if tags else "")


# protein cell colours: 0 gap, 1 match-focal, 2 mismatch, 3 stop(*)
PROT_CMAP = ListedColormap(["#ffffff", "#2c7fb8", "#c7e0ee", "#d7301f"])
# dna cell colours: 0 gap, 1 telomeric, 2 other base
DNA_CMAP = ListedColormap(["#ffffff", "#e6701a", "#c9c9c9"])


DNA_COL = {"telo": "#e6701a", "gap": "#cccccc", "base": "#333333"}


def render_page(pdf, tid, gene, focal, residue, prot_path, dna_path, f5_path, f3_path,
                window, max_rows, motif, block=80, dna_rows=10, max_blocks=16,
                min_flank_cov=0.7, min_flank_ident=0.45, drop_log=None):
    # ---- analyze protein MSA first. An ortholog is trusted only if BOTH flanking
    #      exons (5' AND 3' of the intron) align well -- enough coverage AND enough
    #      identity to the focal. Anything failing on either flank is dropped from ALL
    #      panels (its intron-position / flank / intron-DNA calls are untrustworthy). ----
    prot_rows, p_lo, p_hi, p_c0, focal_res = [], 0, 0, 0, ""
    drop_oids = set()
    if prot_path and os.path.exists(prot_path):
        aln = read_aln(prot_path)
        foc_i = next((i for i, (n, _) in enumerate(aln) if "FOCAL" in n), 0)
        fseq = aln[foc_i][1]
        p_c0 = focal_col(fseq, residue)
        p_lo, p_hi = max(0, p_c0 - window), min(len(fseq), p_c0 + window)
        focal_res = fseq[p_lo:p_hi]

        def flank_quality(seq, cols):
            anchored = [c for c in cols if c < len(fseq) and fseq[c] != "-"]
            if not anchored:
                return 1.0, 1.0
            aligned = [c for c in anchored if c < len(seq) and seq[c] != "-"]
            cov = len(aligned) / len(anchored)
            ident = (sum(1 for c in aligned if seq[c] == fseq[c]) / len(aligned)) if aligned else 0.0
            return cov, ident

        up_cols, dn_cols = range(p_lo, p_c0), range(p_c0, p_hi)
        kept = []
        for i, (nm, seq) in enumerate(aln):
            if i == foc_i:
                continue
            uc, ui = flank_quality(seq, up_cols)
            dc, di = flank_quality(seq, dn_cols)
            # require BOTH flanks: coverage AND identity above thresholds
            if min(uc, dc) < min_flank_cov or min(ui, di) < min_flank_ident:
                oid = nm.split("|")[0]
                drop_oids.add(oid)
                # sidecar: telotron, ortholog_id, up_cov, up_ident, dn_cov, dn_ident,
                #          min_flank_cov_thr, min_flank_ident_thr, reason
                if drop_log is not None:
                    reason = "low_cov" if min(uc, dc) < min_flank_cov else "low_ident"
                    drop_log.append((tid, oid, f"{uc:.3f}", f"{ui:.3f}",
                                     f"{dc:.3f}", f"{di:.3f}",
                                     f"{min_flank_cov:.2f}", f"{min_flank_ident:.2f}",
                                     reason))
            else:
                kept.append((nm, seq))
        kept.sort(key=lambda x: -(float(m.group(1)) if (m := re.search(r"id=([0-9.]+)", x[0])) else 0))
        prot_rows = [aln[foc_i]] + kept[:max_rows - 1]

    # ---- flanking-exon DNA alignments (5' and 3'), poor-upstream orthologs dropped ----
    flank_rows = []   # (label, seq5, seq3); focal first
    if f5_path and os.path.exists(f5_path) and f3_path and os.path.exists(f3_path):
        d5, order = {}, []
        for nm, s in read_aln(f5_path):
            k = nm.split("|")[0]; d5[k] = (nm, s)
            if k not in order:
                order.append(k)
        d3 = {nm.split("|")[0]: (nm, s) for nm, s in read_aln(f3_path)}
        for k in order:
            if k != focal and k in drop_oids:
                continue
            n5, s5 = d5[k]; s3 = d3.get(k, ("", ""))[1]
            flank_rows.append((n5, s5, s3))
        flank_rows = flank_rows[:max_rows]

    # ---- DNA alignment: drop the same poor-upstream orthologs, then size the page ----
    dna_rows_list, dna_W, n_blocks, dna_masks = [], 0, 0, []
    if dna_path and os.path.exists(dna_path):
        aln = read_aln(dna_path)
        foc_i = next((i for i, (n, _) in enumerate(aln) if "TELOTRON" in n), 0)
        others = [a for i, a in enumerate(aln)
                  if i != foc_i and a[0].split("|")[0] not in drop_oids]
        others.sort(key=lambda x: -(float(m.group(1)) if (m := re.search(r"telo=([0-9.]+)", x[0])) else 0))
        dna_rows_list = [aln[foc_i]] + others[:dna_rows - 1]
        dna_W = min(max(len(s) for _, s in dna_rows_list), block * max_blocks)
        dna_masks = [telomeric_mask(s, motif) for _, s in dna_rows_list]
        n_blocks = (dna_W + block - 1) // block

    # ---- fixed character-cell geometry: uniform, legible physical text size ----
    CW = 0.12                      # inches per alignment column (all panels)
    RH, RHI = 0.205, 0.175         # row height: protein/flank, intron
    FS, FSI, LF, TF, SF = 9.0, 7.5, 8.5, 12, 14   # char, intron-char, label, title, suptitle
    LBLW = 20                      # in-axes columns reserved at left for row labels
    FLANK_SHOW, SEP = 36, 4        # junction-proximal bp of each flank shown; gap cols

    n_prot, n_fl, n_dna = len(prot_rows), len(flank_rows), len(dna_rows_list)
    prot_cols = (p_hi - p_lo) if prot_rows else 0
    fullW5 = max((len(s) for _, s, _ in flank_rows), default=0)
    fullW3 = max((len(s) for _, _, s in flank_rows), default=0)
    show5, show3 = min(fullW5, FLANK_SHOW), min(fullW3, FLANK_SHOW)
    flank_cols = (show5 + SEP + show3) if flank_rows else 0

    prot_h = n_prot * RH + 0.55 if prot_rows else 0.3
    flank_h = n_fl * RH + 0.8 if flank_rows else 0.3
    dna_h = n_blocks * (n_dna * RHI + 0.42) + 0.35 if dna_rows_list else 0.3
    GAP, TOP, BOT, LEFTIN, RIGHTIN = 0.5, 1.0, 0.2, 0.15, 0.35
    pw = {"p": (LBLW + prot_cols) * CW, "f": (LBLW + flank_cols) * CW, "d": (LBLW + block) * CW}
    Wf = LEFTIN + max(pw.values()) + RIGHTIN
    Hf = TOP + BOT + prot_h + flank_h + dna_h + 2 * GAP
    fig = plt.figure(figsize=(Wf, Hf))
    fig.suptitle(f"{gene}    {tid}    intron @ residue {residue}", fontsize=SF, x=0.012, ha="left",
                 y=1 - 0.3 / Hf)
    cur = Hf - TOP
    axp = fig.add_axes([LEFTIN / Wf, (cur - prot_h) / Hf, pw["p"] / Wf, prot_h / Hf]); cur -= prot_h + GAP
    axf = fig.add_axes([LEFTIN / Wf, (cur - flank_h) / Hf, pw["f"] / Wf, flank_h / Hf]); cur -= flank_h + GAP
    axd = fig.add_axes([LEFTIN / Wf, (cur - dna_h) / Hf, pw["d"] / Wf, dna_h / Hf])

    # ---- A: protein MSA window around the intron ----
    if prot_rows:
        lo, hi, c0, ncol = p_lo, p_hi, p_c0, p_hi - p_lo
        mat = np.zeros((n_prot, ncol))
        for ri, (nm, seq) in enumerate(prot_rows):
            for j, col in enumerate(range(lo, hi)):
                ch = seq[col] if col < len(seq) else "-"
                mat[ri, j] = (0 if ch == "-" else 3 if ch == "*"
                              else 1 if (ch == focal_res[j] and focal_res[j] != "-") else 2)
        axp.imshow(mat, aspect="auto", cmap=PROT_CMAP, vmin=0, vmax=3,
                   interpolation="nearest", extent=[0, ncol, n_prot, 0])
        for ri, (nm, seq) in enumerate(prot_rows):
            for j, col in enumerate(range(lo, hi)):
                ch = seq[col] if col < len(seq) else "-"
                if ch != "-":
                    axp.text(j + 0.5, ri + 0.5, ch, ha="center", va="center", family="monospace",
                             fontsize=FS, color="white" if mat[ri, j] in (1, 3) else "black")
            axp.text(-1, ri + 0.5, short(nm), ha="right", va="center", fontsize=LF)
        axp.axvline((c0 - lo) + 1, color="#111111", lw=2.4)
        for col in range(lo, hi):
            if col % 10 == 0:
                axp.text(col - lo + 0.5, n_prot + 0.5, str(col), ha="center", va="top",
                         fontsize=6.5, color="#888")
        axp.set_xlim(-LBLW, ncol); axp.set_ylim(n_prot + 1.05, -0.1); axp.axis("off")
        ttl = "A  flanking-exon protein alignment   (black line = intron)"
        if drop_oids:
            ttl += f"   —   {len(drop_oids)} ortholog(s) dropped (poor flank alignment)"
        axp.set_title(ttl, fontsize=TF, loc="left")
    else:
        axp.text(0.5, 0.5, "no protein MSA", ha="center", va="center"); axp.axis("off")

    # ---- B: flanking-exon DNA alignment (5' | intron | 3'), junction-proximal ----
    if flank_rows:
        foc5, foc3 = flank_rows[0][1], flank_rows[0][2]

        def draw(seqs, focal_seq, cols, x0):
            for ri, s in enumerate(seqs):
                for k, col in enumerate(cols):
                    ch = s[col] if col < len(s) else "-"
                    if ch == "-":
                        cc, txt = DNA_COL["gap"], "-"
                    elif col < len(focal_seq) and ch == focal_seq[col]:
                        cc, txt = "#1a7a3a", ch
                    else:
                        cc, txt = "#999999", ch
                    axf.text(x0 + k, ri, txt, ha="center", va="center", family="monospace",
                             fontsize=FS, color=cc)
        draw([s for _, s, _ in flank_rows], foc5, range(fullW5 - show5, fullW5), 0)
        draw([s for _, _, s in flank_rows], foc3, range(0, show3), show5 + SEP)
        axf.axvspan(show5 - 0.5, show5 + SEP - 0.5, color="#ede7f6")
        axf.text(show5 + SEP / 2 - 0.5, -1.0, "intron", ha="center", va="bottom", fontsize=7.5, color="#5e35b1")
        for ri, (lab, _, _) in enumerate(flank_rows):
            axf.text(-1, ri, short(lab), ha="right", va="center", fontsize=LF)
        axf.text(0, n_fl - 0.25, "5' flank exon →| donor", ha="left", va="top", fontsize=7.5)
        axf.text(show5 + SEP, n_fl - 0.25, "acceptor |→ 3' flank exon", ha="left", va="top", fontsize=7.5)
        axf.set_xlim(-LBLW, show5 + SEP + show3); axf.set_ylim(n_fl + 0.5, -1.4); axf.axis("off")
        axf.set_title("B  flanking-exon DNA alignment   (green = matches focal)", fontsize=TF, loc="left")
    else:
        axf.text(0.5, 0.5, "no flank alignment", ha="center", va="center"); axf.axis("off")

    # ---- C: telotron vs orthologous-intron DNA, wrapped characters ----
    if dna_rows_list:
        labels = [short(n) for n, _ in dna_rows_list]
        for b in range(n_blocks):
            y0 = b * (n_dna + 1.5)
            cstart, cend = b * block, min(b * block + block, dna_W)
            for ri, (nm, seq) in enumerate(dna_rows_list):
                y = y0 + ri; tm = dna_masks[ri]
                for col in range(cstart, cend):
                    ch = seq[col] if col < len(seq) else "-"
                    if ch == "-":
                        cc, txt = DNA_COL["gap"], "-"
                    elif col < len(tm) and tm[col]:
                        cc, txt = DNA_COL["telo"], ch
                    else:
                        cc, txt = DNA_COL["base"], ch
                    axd.text(col - cstart, y, txt, ha="center", va="center", family="monospace",
                             fontsize=FSI, color=cc)
                axd.text(-1, y, labels[ri], ha="right", va="center", fontsize=LF)
            axd.text(-1, y0 - 0.6, f"{cstart + 1}", ha="right", va="center", fontsize=6.5, color="#888")
        axd.set_xlim(-LBLW, block); axd.set_ylim(n_blocks * (n_dna + 1.5), -1.1); axd.axis("off")
        trunc = "   (truncated)" if (max(len(s) for _, s in dna_rows_list) > dna_W) else ""
        axd.set_title(f"C  telotron vs orthologous-intron DNA   (orange = telomeric repeat){trunc}",
                      fontsize=TF, loc="left")
    else:
        axd.text(0.5, 0.5, "no orthologous intron found (CREATE / no-align)",
                 ha="center", va="center"); axd.axis("off")

    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ortho-dir", required=True, help="telotron_ortholog_align.py output dir")
    ap.add_argument("--out", required=True, help="output multipage PDF")
    ap.add_argument("--window", type=int, default=28, help="+/- residues of protein flank to show")
    ap.add_argument("--max-rows", type=int, default=12, help="max sequences per panel")
    ap.add_argument("--motif", default="TTTAGGG")
    ap.add_argument("--max-loci", type=int, default=0, help="0 = all")
    ap.add_argument("--min-flank-cov", type=float, default=0.7,
                    help="drop orthologs unless BOTH flanks align with >= this coverage")
    ap.add_argument("--min-flank-ident", type=float, default=0.45,
                    help="drop orthologs unless BOTH flanks reach >= this identity to focal")
    ap.add_argument("--png-dir", default=None)
    ap.add_argument("--drop-log", default=None,
                    help="optional TSV: dropped orthologs (telotron, ortholog_id, "
                         "up_cov, up_ident, dn_cov, dn_ident, thr_cov, thr_ident, reason); "
                         "defaults to <out>.dropped_orthologs.tsv alongside --out")
    args = ap.parse_args()

    # one record per telotron from ortholog_map.tsv (tid -> focal, gene, protein, residue)
    loci = {}
    mp = os.path.join(args.ortho_dir, "ortholog_map.tsv")
    for line in open(mp):
        f = line.rstrip("\n").split("\t")
        if f[0] == "telotron" or len(f) < 5:
            continue
        loci.setdefault(f[0], dict(focal=f[1], gene=f[2], protein=f[3], residue=int(f[4])))
    order = sorted(loci.items(), key=lambda kv: (kv[1]["focal"], kv[1]["gene"]))
    if args.max_loci:
        order = order[:args.max_loci]
    if args.png_dir:
        os.makedirs(args.png_dir, exist_ok=True)

    pmsa = os.path.join(args.ortho_dir, "protein_msa")
    idna = os.path.join(args.ortho_dir, "intron_dna")
    fdna = os.path.join(args.ortho_dir, "flank_dna")
    drop_log = []
    n = 0
    with PdfPages(args.out) as pdf:
        for tid, d in order:
            prot_path = os.path.join(pmsa, f"{d['focal']}__{d['protein']}.aln.fasta")
            safe = tid.replace("|", "__").replace(":", "_")
            dna_path = os.path.join(idna, f"{safe}.aln.fasta")
            f5_path = os.path.join(fdna, f"{safe}.5p.aln.fasta")
            f3_path = os.path.join(fdna, f"{safe}.3p.aln.fasta")
            if not os.path.exists(prot_path) and not os.path.exists(dna_path):
                continue
            render_page(pdf, tid, d["gene"], d["focal"], d["residue"],
                        prot_path, dna_path, f5_path, f3_path, args.window, args.max_rows, args.motif,
                        min_flank_cov=args.min_flank_cov, min_flank_ident=args.min_flank_ident,
                        drop_log=drop_log)
            n += 1
            if n % 50 == 0:
                print(f"  {n} pages...", file=sys.stderr)
    drop_path = args.drop_log or (args.out + ".dropped_orthologs.tsv")
    with open(drop_path, "w") as fh:
        fh.write("telotron\tortholog_id\tup_cov\tup_ident\tdn_cov\tdn_ident\t"
                 "thr_cov\tthr_ident\treason\n")
        for row in drop_log:
            fh.write("\t".join(row) + "\n")
    print(f"wrote {args.out} ({n} locus pages); "
          f"dropped {len(drop_log)} ortholog rows -> {drop_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
